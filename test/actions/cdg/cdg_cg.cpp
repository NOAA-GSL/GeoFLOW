//==================================================================================
// Module       : gtest_cg.cpp
// Date         : 3/16/20 (DLR)
// Description  : GeoFLOW test of conjigate gradient operator,
//                solving
//                    Nabla^2 u = -f,
//                In a 2D box, where f =  6xy(1-y) - 2x^3
//                on 0 <= x,y <= 1, and u(x,y=0)=u(x,y=1)=0;
//                u(x=0,y)=0; u(x=1,y) = y(1-y).
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

#include <unistd.h>

#include <cstdio>
#include <iostream>

#include "gcg.hpp"
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "ggrid_box.hpp"
#include "ggrid_factory.hpp"
#include "ghelmholtz.hpp"
#include "gio.hpp"
#include "gio_observer.hpp"
#include "gllbasis.hpp"
#include "gmass.hpp"
#include "gmorton_keygen.hpp"
#include "gtypes.h"
#include "pdeint/io_base.hpp"
#include "pdeint/null_observer.hpp"
#include "pdeint/observer_base.hpp"
#include "pdeint/observer_factory.hpp"
#include "tbox/error_handler.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/property_tree.hpp"

using namespace geoflow::pdeint;
using namespace geoflow::tbox;
using namespace std;

struct stStateInfo {
    GINT sttype = 1;    // state type index (grid=0 or state=1)
    GINT gtype = 0;     // check src/cdg/include/gtypes.h
    GSIZET index = 0;   // time index
    GSIZET nelems = 0;  // num elems
    GSIZET cycle = 0;   // continuous time cycle
    GFTYPE time = 0.0;  // state time
    std::vector<GString>
        svars;  // names of state members
    GTVector<GStateCompType>
        icomptype;  // encoding of state component types
    GTMatrix<GINT>
        porder;    // if ivers=0, is 1 X GDIM; else nelems X GDIM;
    GString idir;  // input directory
    GString odir;  // output directory
};

struct TypePack {
    using State = GTVector<GTVector<GFTYPE> *>;
    using StateComp = GTVector<GFTYPE>;
    using StateInfo = GStateInfo;
    using Grid = GGrid<TypePack>;
    using GridBox = GGridBox<TypePack>;
    using GridIcos = GGridIcos<TypePack>;
    using Mass = GMass<TypePack>;
    using Ftype = GFTYPE;
    using Derivative = State;
    using Time = Ftype;
    using CompDesc = GTVector<GStateCompType>;
    using Jacobian = State;
    using Size = GSIZET;
    using EqnBase = EquationBase<TypePack>;       // Equation Base type
    using EqnBasePtr = std::shared_ptr<EqnBase>;  // Equation Base ptr
    using IBdyVol = GTVector<GSIZET>;
    using TBdyVol = GTVector<GBdyType>;
    using Operator = GHelmholtz<TypePack>;
    using GElemList = GTVector<GElem_base *>;
    using Preconditioner = GHelmholtz<TypePack>;
    using ConnectivityOp = GGFX<Ftype>;
    using FilterBasePtr = std::shared_ptr<FilterBase<TypePack>>;
    using FilterList = std::vector<FilterBasePtr>;
};
struct CGTypes {
    using Operator = class GHelmholtz<TypePack>;
    using Preconditioner = GLinOpBase<TypePack>;
    using State = typename TypePack::State;
    using StateComp = typename TypePack::StateComp;
    using Grid = GGrid<TypePack>;
    using Ftype = typename TypePack::Ftype;
    using ConnectivityOp = GGFX<GFTYPE>;
};

using Types = TypePack;                         // Define types used
using EqnBase = EquationBase<Types>;            // Equation Base type
using EqnBasePtr = std::shared_ptr<EqnBase>;    // Equation Base ptr
using IOBaseType = IOBase<Types>;               // IO Base type
using IOBasePtr = std::shared_ptr<IOBaseType>;  // IO Base ptr
using Grid = TypePack::Grid;
using Ftype = TypePack::Ftype;

Grid *grid_ = NULLPTR;
GC_COMM comm_ = GC_COMM_WORLD;  // communicator

GINT szMatCache_ = _G_MAT_CACHE_SIZE;
GINT szVecCache_ = _G_VEC_CACHE_SIZE;

void init_ggfx(PropertyTree &ptree, Grid &grid, GGFX<Ftype> &ggfx);

int main(int argc, char **argv) {
    GString serr = "main: ";
    GINT errcode, gerrcode, iret;
    IOBasePtr pIO;  // ptr to IOBase operator
    GTVector<GTVector<Ftype>> *xnodes;
    LinSolverBase<CGTypes>::Traits cgtraits;
    GString snorm;
    std::ifstream itst;
    std::ofstream ios;

    typename ObserverBase<Types>::Traits
        binobstraits;

    if (argc > 1) {
        cout << "No arguments accepted" << endl;
        exit(1);
    }
    errcode = 0;

    // Initialize comm:
    GComm::InitComm(&argc, &argv);
    mpixx::environment env(argc, argv);  // init GeoFLOW comm
    mpixx::communicator world;
    GlobalManager::initialize(argc, argv);
    GlobalManager::startup();

    GINT myrank = GComm::WorldRank();
    GINT nprocs = GComm::WorldSize();

    // Get minimal property tree:
    PropertyTree gtree, ptree;
    GString sgrid;
    GTPoint<Ftype> P0, P1, dP;
    std::vector<GINT> pstd(GDIM);

    ptree.load_file("cg_input.jsn");
    sgrid = ptree.getValue<GString>("grid_type");
    pstd = ptree.getArray<GINT>("exp_order");
    gtree = ptree.getPropertyTree(sgrid);

    std::vector<Ftype> xyz0 = gtree.getArray<Ftype>("xyz0");
    std::vector<Ftype> dxyz = gtree.getArray<Ftype>("delxyz");
    P0 = xyz0;
    P1 = dxyz;
    P1 += P0;

    assert(P0.x1 == 0.0 && P1.x1 == 1.0 && P0.x2 == 0.0 && P1.x2 == 1.0);

    // Create basis:
    GTVector<GNBasis<GCTYPE, Ftype> *> gbasis(GDIM);
    for (GSIZET k = 0; k < GDIM; k++) {
        gbasis[k] = new GLLBasis<GCTYPE, Ftype>(pstd[k]);
    }

    EH_MESSAGE("main: Create grid...");

    ObserverFactory<Types>::get_traits(ptree, "gio_observer", binobstraits);
    grid_ = GGridFactory<Types>::build(ptree, gbasis, pIO, binobstraits, comm_);

    EH_MESSAGE("main: Initialize GGFX operator...");

    // Initialize gather/scatter operator:
    GGFX<Ftype> ggfx;
    init_ggfx(ptree, *grid_, ggfx);

    xnodes = &(grid_->xNodes());

    GTVector<Ftype> f(grid_->ndof());
    GTVector<Ftype> u(grid_->ndof());
    GTVector<Ftype> ua(grid_->ndof());
    GTVector<Ftype> ub(grid_->ndof());
    GTVector<Ftype> *mass_local = grid_->massop().data();
    GTVector<GTVector<Ftype> *> utmp(10);

    for (auto j = 0; j < utmp.size(); j++) utmp[j] = new GTVector<Ftype>(grid_->ndof());

    cgtraits.maxit = gtree.getValue<GDOUBLE>("maxit");
    cgtraits.tol = gtree.getValue<GDOUBLE>("tol");
    snorm = gtree.getValue<GString>("norm_type");
    cgtraits.normtype = LinSolverBase<CGTypes>::str2normtype(snorm);

    EH_MESSAGE("main: Initialize CG operator...");

    // Initialize GCG operator:
    GSIZET niter;
    Ftype eps = 100.0 * std::numeric_limits<Ftype>::epsilon();
    Ftype err, resmin, resmax, x, y, z;
    GTVector<Ftype> *resvec;

    EH_MESSAGE("main: Create Lap op...");
    GHelmholtz<Types> L(*grid_);  // Laplacian operator

    EH_MESSAGE("main: Create CG op...");

    GCG<CGTypes> cg(cgtraits, *grid_, ggfx, utmp);

    EH_MESSAGE("main: Check for SPD...");

    //L.use_metric(FALSE);

    // Check if operator is SPD:
    GTMatrix<Ftype> Hmat(f.size(), f.size());
    f = 0.0;
    for (auto j = 0; j < f.size(); j++) {
        f[j] = 1.0;
        L.opVec_prod(f, utmp, u);
        for (auto i = 0; i < f.size(); i++) Hmat(i, j) = u[i];
        f[j] = 0.0;
    }

//cout << "Hmat=" << Hmat << endl;
#if 1
    GSIZET nbad = 0, ngbad;
    for (auto j = 0; j < f.size(); j++) {
        for (auto i = j; i < f.size(); i++) {
            if (!FUZZYEQ(Hmat(i, j), Hmat(j, i), eps)) {
                cout << "main: (" << i << "," << j << "): H=" << Hmat(i, j) << " H^T=" << Hmat(j, i) << endl;
                nbad++;
            }
        }
    }

    // Accumulate 'bad' entry number:
    GComm::Allreduce(&nbad, &ngbad, 1, T2GCDatatype<GSIZET>(), GC_OP_MAX, comm_);
    if (ngbad == 0) {
        cout << "main: .................................. operator is SPD!" << endl;
    } else {
        cout << "main: .................................. operator NOT SPD!" << endl;
        cout << "main: .................................. ngbad=" << nbad << "; size=" << f.size() * (f.size() + 1) / 2 << endl;
        errcode = 1;
        goto prerror;
    }
#endif

    EH_MESSAGE("main: Set RHS...");

    // Generate smooth RHS, bdyy vectors:
    //    f =  6xy(1-y) - 2x^3

    ub = 0.0;
    for (auto j = 0; j < grid_->ndof(); j++) {
        x = (*xnodes)[0][j];
        y = (*xnodes)[1][j];
        if (GDIM > 2) z = (*xnodes)[2][j];
        f[j] = 6.0 * x * y * (1.0 - y) - 2.0 * pow(x, 3.0);    // RHS
        ua[j] = y * (1.0 - y) * pow(x, 3.0);                   // analytic solution
        if (FUZZYEQ(P0.x2, y, eps) || FUZZYEQ(P1.x2, y, eps))  // N & S bdy
            ub[j] = 0.0;
        if (FUZZYEQ(P0.x1, x, eps))  // W bdy
            ub[j] = 0.0;
        if (FUZZYEQ(P1.x1, x, eps))  // E bdy
            ub[j] = y * (1.0 - y);
    }

    // Multiply f by local mass matrix:
    f.pointProd(*mass_local);  // M_L f_L
    f *= -1.0;                 // -f required by discretization

    EH_MESSAGE("main: Solve linear system...");

    // Invert mass operator, solve L u = f for u,
    // where L is Laplacian operator:
    u = 0.0;  // initial guess
    iret = cg.solve(L, f, ub, u);

    niter = cg.get_iteration_count();
    resmin = cg.get_resid_min();
    resmax = cg.get_resid_max();
    resvec = &cg.get_residuals();
    EH_MESSAGE("main: Compute errors...");

    *utmp[0] = u - ua;
    utmp[0]->rpow(2);
    err = grid_->integrate(*utmp[0], *utmp[1]) * grid_->ivolume();

    EH_MESSAGE("main: Write to file...");

    itst.open("cg_err.txt");
    ios.open("cg_err.txt", std::ios_base::app);

    if (itst.peek() == std::ofstream::traits_type::eof()) {
        ios << "#elems"
            << "  "
            << "#dof"
            << "  ";
        for (auto j = 0; j < GDIM; j++) ios << "p" << j + 1 << "  ";
        ios << "L2_err"
            << "  "
            << "niter"
            << "  "
            << "resid_min"
            << "  "
            << "resid_max" << std::endl;
    }
    itst.close();

    ios << grid_->ngelems() << "  "
        << grid_->ndof() << "  ";
    for (auto j = 0; j < GDIM; j++)
        ios << pstd[j] << "  ";
    ios << err << "  ";
    ios << niter << "  ";
    ios << resmin << "  ";
    ios << resmax << std::endl;

    ios.close();

    errcode = err < 1e-12 ? 0 : 2;
    if (errcode != 0 || iret != GCG<CGTypes>::GCGERR_NONE) {
        cout << serr << " Error: err=" << err << " code=" << errcode << endl;
        cout << serr << " residuals=" << *resvec << endl;
    }
    assert(iret == GCG<CGTypes>::GCGERR_NONE && "Solve failure");

prerror:
    // Accumulate error codes:
    GComm::Allreduce(&errcode, &gerrcode, 1, T2GCDatatype<GINT>(), GC_OP_MAX, comm_);

    if (gerrcode != 0) {
        cout << serr << " Error: errcode=" << gerrcode << endl;
    } else {
        cout << serr << " Success!" << endl;
    }

    GComm::TermComm();

    return (gerrcode);

}  // end, main

//==================================================================================
// Module       : gtools.cpp
// Date         : 5/5/19 (DLR)
// Description  : GeoFLOW initialization tools
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From :
//==================================================================================
//#include "gtools.h"

//**********************************************************************************
//**********************************************************************************
// METHOD: init_ggfx
// DESC  : Initialize gather/scatter operator
// ARGS  : ptree   : main property tree
//         grid    : Grid object
//         ggfx    : gather/scatter op, GGFX
//**********************************************************************************
//**********************************************************************************
void init_ggfx(PropertyTree &ptree, Grid &grid, GGFX<Ftype> &ggfx) {
    const auto ndof = grid_->ndof();
    const auto nxyz = grid.xNodes().size();
    ASSERT(nxyz <= GGFX<Ftype>::NDIM);
    std::vector<std::array<Ftype, GGFX<Ftype>::NDIM>> xyz(ndof);
    for (std::size_t i = 0; i < ndof; i++) {
        for (std::size_t d = 0; d < nxyz; d++) {
            xyz[i][d] = grid.xNodes()[d][i];
        }
        for (std::size_t d = nxyz; d < GGFX<Ftype>::NDIM; d++) {
            xyz[i][d] = 0;
        }
    }

    // Get max duplicate points
    std::size_t maxdups = 0;
    switch (GDIM) {
        case 1:
            maxdups = 2;
            break;
        case 2:
            maxdups = 6;
            break;
        case 3:
            maxdups = 12;
            break;
    }

    // Create GGFX
    ggfx.init(maxdups, 0.25 * grid.minnodedist(), xyz);

}  // end method init_ggfx
