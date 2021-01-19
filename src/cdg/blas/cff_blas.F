c======================================================================================
c Name         : cff_blas
c Date         : 1/1/18 (DLR)
c Copyright    : Copyright 2018. Colorado State University. All rights reserved
c Description  : Cache-friendly optimized Fortran linear algebra routines
c======================================================================================


c======================================================================================
c Quad-precision versions
c======================================================================================


c************************************************************************************
c************************************************************************************
c METHOD     : qmxm
c DESCRIPTION: mxm driver routine that that calls methods explicitly unrolled
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm(C,A,nai,naj,B,nbi,nbj,isz) bind(C,name="qmxm")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      integer     (C_INT), intent (in) :: isz
      integer  (C_SIZE_T), intent (in) :: nai, naj, nbi, nbj
      real(C_LONG_DOUBLE), intent (in) :: A(nai,naj), B(nbi,nbj)
      real(C_LONG_DOUBLE), intent(out) :: C(nai,nbj)
      real(C_LONG_DOUBLE)              :: alpha, beta

!!$acc routine vector qmxm1
!!$acc routine vector qmxm2
!!$acc routine vector qmxm3
!!$acc routine vector qmxm4
!!$acc routine vector qmxm5
!!$acc routine vector qmxm6
!!$acc routine vector qmxm7
!!$acc routine vector qmxm8
!!$acc routine vector qmxm9
!!$acc routine vector qmxm10
!!$acc routine vector qmxm11
!!$acc routine vector qmxm12
!!$acc routine vector qmxm13
!!$acc routine vector qmxm14
!!$acc routine vector qmxm15
!!$acc routine vector qmxm16
!!$acc routine vector qmxm17
!!$acc routine vector qmxm18
!!$acc routine vector qmxm19
!!$acc routine vector qmxm20
!!$acc routine vector qmxm21
!!$acc routine vector qmxm22
!!$acc routine vector qmxm23
!!$acc routine vector qmxm24
!!!$acc routine vector qmxmcf
      if      ( nbi .eq.  1 ) then
        call qmxm1 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  2 ) then
        call qmxm2 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  3 ) then
        call qmxm3 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  4 ) then
        call qmxm4 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  5 ) then
        call qmxm5 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  6 ) then
        call qmxm6 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  7 ) then
        call qmxm7 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  8 ) then
        call qmxm8 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  9 ) then
        call qmxm9 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 10 ) then
        call qmxm10(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 11 ) then
        call qmxm11(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 12 ) then
        call qmxm12(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 13 ) then
        call qmxm13(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 14 ) then
        call qmxm14(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 15 ) then
        call qmxm15(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 16 ) then
        call qmxm16(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 17 ) then
        call qmxm17(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 18 ) then
        call qmxm18(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 19 ) then
        call qmxm19(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 20 ) then
        call qmxm20(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 21 ) then
        call qmxm21(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 22 ) then
        call qmxm22(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 23 ) then
        call qmxm23(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 24 ) then
        call qmxm24(A, nai, B, nbj, C, nbj)
      else
        call qmxmcf(C, A, nai, naj, B, nbi, nbj, isz )
#if 0
        alpha = 1.0
        beta  = 0.0
        call dgemm("N","N", nai, nbj, alpha, A, nai, 
     1             B, nbi, beta, C, nai)
#endif
      endif


      return
c
c  end of subroutine qmxm
c
      end subroutine qmxm


c************************************************************************************
c************************************************************************************
c METHOD     : qmxmcf
c DESCRIPTION: Pure Fortran cache-friendly mxm routine that does not unroll by hand,
c              but instead uses cache blocking for mat-mat multiply
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxmcf(C,A,nai,naj,B,nbi,nbj,isz)
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer  (C_SIZE_T), intent   (in) :: nai, naj, nbi, nbj
      integer     (C_INT), intent   (in) :: isz
      real(C_LONG_DOUBLE), intent   (in) :: A(nai,naj), B(nbi,nbj)
      real(C_LONG_DOUBLE), intent  (out) :: C(nai,nbj)
      real(C_LONG_DOUBLE)                :: zero
      integer  (C_SIZE_T)                :: i, ii, j, jj, k, kk

c  Initialize C to 0:
      !$acc parallel

!     zero = 0.0
!     call matset(C, nai, nbj, zero, isz)
!     'matset' seems to be required to be provided here 
!     explicitly as ACC isn't finding its device code:
      do ii = 1, nai, isz
        do jj = 1, nbj, isz

          do i = ii, min(nai,ii+isz-1)
            do j = jj, min(nbj,jj+isz-1)
                C(i,j) = 0.0
            enddo
          enddo

        enddo
      enddo


c  Do matrix-matrix multiply using cache-blocking:
      do ii = 1, nai, isz
        do jj = 1, nbj, isz
          do kk = 1, nbi, isz
  
            do i = ii, min(nai,ii+isz-1)
              do j = jj, min(nbj,jj+isz-1)
                do k = kk, min(nbi,kk+isz-1)
                  C(i,j) = C(i,j) + A(i,k)*B(k,j)
                enddo
              enddo
            enddo
        
          enddo
        enddo
      enddo

      !$acc end parallel

      return
c
c  end of subroutine qmxmcf
c
      end subroutine qmxmcf

c************************************************************************************
c************************************************************************************
c METHOD     : qmxv
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxv(y, A, x, nx, ny, isz)  bind(C,name="qmxv")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent  (in)     :: isz
      integer  (C_SIZE_T), intent  (in)     :: ny, nx
      real(C_LONG_DOUBLE), intent  (in)     :: A(ny,*), x(*)
      real(C_LONG_DOUBLE), intent (out)     :: y(*)
      integer  (C_SIZE_T)                   :: i, ii, j, jj

      do ii = 1, ny, isz
        do i = ii, min(ny,ii+isz-1)
           y(i) = 0.0
        enddo
      enddo


c  Do matrix-vector multiply:
      do ii = 1, ny, isz
        do jj = 1, nx, isz

          do i = ii, min(ny,ii+isz-1)
            do j = jj, min(nx,jj+isz-1)
              y(i) = y(i) + A(i,j)*x(j)
            enddo
          enddo

        enddo
      enddo

      return
c
c  end of subroutine qmxv
c
      end subroutine qmxv

c************************************************************************************
c************************************************************************************
c METHOD     : qmxDm  
c DESCRIPTION: computes C = A * Diag(B). Dimensions of C assumed to be 
c              those of A
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxDm(C,A,nai,naj,b,nb,isz) bind(C,name="qmxDm")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent (in)  ::  isz
      integer  (C_SIZE_T), intent (in)  ::  nai, naj, nb
      real(C_LONG_DOUBLE), intent (in)  ::  A(nai,naj), b(*)
      real(C_LONG_DOUBLE), intent(out)  ::  C(nai,naj)

      integer(C_SIZE_T)               ::  i, ii, j, jj

c  Do matrix-Diag(matrix) multiply:
      do ii = 1, nai, isz 
        do jj = 1, naj, isz 

          do i = ii, min(nai,ii+isz-1)
            do j = jj, min(naj,jj+isz-1)
              C(i,j) = A(i,j)*b(j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine qmxDm
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : qDmxm  
c DESCRIPTION: computes C = Diag(A)*B. Dimensions of C assumed to be 
c              those of B
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qDmxm(C, a, na, B, nbi, nbj, isz) bind(C,name="qDmxm")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent (in) ::  isz
      integer  (C_SIZE_T), intent (in) ::  nbi, nbj, na
      real(C_LONG_DOUBLE), intent (in) ::  B(nbi,nbj), a(na)
      real(C_LONG_DOUBLE), intent(out) ::  C(nbi,nbj)
      integer  (C_SIZE_T)              ::  i, ii, j, jj

c  Do Diag(matrix)-matrix multiply:
      do ii = 1, nbi, isz 
        do jj = 1, nbj, isz 

          do i = ii, min(nbi,ii+isz-1)
            do j = jj, min(nbj,jj+isz-1)
              C(i,j) = a(i)*B(i,j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine qDmxm
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : qaApbB  
c DESCRIPTION: computes C = a A + b B. Dimensions of A, B, C assumed to be equal
c              (nxm), and a, and b are double constants
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qaApbB(C,A,B,n,m,ca,cb,isz) bind(C,name="qaApbB")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent (in) :: isz
      integer  (C_SIZE_T), intent (in) :: n, m
      real(C_LONG_DOUBLE), intent (in) :: A(n,m), B(n,m), ca, cb
      real(C_LONG_DOUBLE), intent(out) :: C(n,m)
      integer  (C_SIZE_T)              :: i, ii, j, jj

c  Do sum:
      do ii = 1, m, isz 
        do jj = 1, n, isz 

          do i = ii, min(m,ii+isz-1)
            do j = jj, min(n,jj+isz-1)
              C(i,j) = ca*A(i,j) + cb*B(i,j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine qaApbB
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : qzaxpby
c DESCRIPTION: z = a*x + b*y
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qzaxpby(z,x,a,y,b,nxy,isz) bind(C,name="qzaxpby")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent (in) ::  isz
      integer  (C_SIZE_T), intent (in) ::  nxy
      real(C_LONG_DOUBLE), intent (in) :: x(nxy), y(nxy), a, b
      real(C_LONG_DOUBLE), intent(out) :: z(nxy)
      real(C_LONG_DOUBLE)              :: prod
      integer  (C_SIZE_T)              :: i, ii
c
c  Do vector sums:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = b*y(i) 
          z(i) = a*x(i) + prod
        enddo
      enddo
#endif

      do i = 1, nxy
        z(i) = a*x(i) + b*y(i)
      enddo

      return
c
c  end of subroutine qzaxpby
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qxaxpby
c DESCRIPTION: x(i) = a*x(i)+b*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qxaxpby(x,a,y,b,nxy,isz) bind(C,name="qxaxpby")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent   (in) :: isz
      integer  (C_SIZE_T), intent   (in) :: nxy
      real(C_LONG_DOUBLE), intent   (in) :: y(*), a, b
      real(C_LONG_DOUBLE), intent(inout) :: x(*)
      real(C_LONG_DOUBLE)                :: prod
      integer  (C_SIZE_T)                :: i, ii

c  Do vector sums:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = b*y(i) + 0.0
          x(i) = a*x(i) + prod
        enddo
      enddo
#endif
      do i = 1, nxy
        x(i) = a*x(i) + b*y(i)
      enddo

      return
c
c  end of subroutine qxaxpby
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qzvxvpt
c DESCRIPTION: z(i) = x(i)*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qzvxvpt(z,x,y,nxy,isz) bind(C,name="qzvxvpt")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent (in) :: isz
      integer  (C_SIZE_T), intent (in) :: nxy
      real(C_LONG_DOUBLE), intent (in) :: x(*), y(*)
      real(C_LONG_DOUBLE), intent(out) :: z(*)
      integer  (C_SIZE_T)              :: i, ii

c  Do point product:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          z(i) = x(i) * y(i) 
        enddo
      enddo
#endif

      do i = 1, nxy
        z(i) = x(i) * y(i) 
      enddo

      return
c
c  end of subroutine qzvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qvvxvpt
c DESCRIPTION: x(i) = x(i)*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qvvxvpt(x,y,nxy,isz) bind(C,name="qvvxvpt")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent   (in) :: isz
      integer  (C_SIZE_T), intent   (in) :: nxy
      real(C_LONG_DOUBLE), intent   (in) :: y(*)
      real(C_LONG_DOUBLE), intent(inout) :: x(*)
      integer  (C_SIZE_T)                :: i, ii

c  Do point product:
#if 0
      do ii = 1, nxy, isz

        do i = ii, min(nxy,ii+isz-1)
          x(i) = x(i) * y(i) 
        enddo

      enddo
#endif

      do i = 1, nxy
        x(i) = x(i) * y(i) 
      enddo

      return
c
c  end of subroutine qvvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qvvxvptpv
c DESCRIPTION: x(i) = x(i)*y(i) + cz*z(i)
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qvvxvptpv(x,y,z,cz,nxy,isz) bind(C,name="qvvxvptpv")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent   (in) :: isz
      integer  (C_SIZE_T), intent   (in) :: nxy
      real(C_LONG_DOUBLE), intent   (in) :: y(*), z(*), cz
      real(C_LONG_DOUBLE), intent(inout) :: x(*)
      real(C_LONG_DOUBLE)                :: prod
      integer  (C_SIZE_T)                :: i, ii

c  Do point product & sum:
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = cz*z(i) 
          x(i) = x(i) * y(i) + prod
        enddo
      enddo

      return
c
c  end of subroutine qvvxvptpv
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qvpvxvpt
c DESCRIPTION: z(i) = x(i)*y(i) + cz*z(i)
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qvpvxvpt(x,y,z,cz,nxy,isz) bind(C,name="qvpvxvpt")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent   (in) :: isz
      integer  (C_SIZE_T), intent   (in) :: nxy
      real(C_LONG_DOUBLE), intent   (in) :: y(*), x(*), cz
      real(C_LONG_DOUBLE), intent(inout) :: z(*)
      real(C_LONG_DOUBLE)                :: prod
      integer  (C_SIZE_T)                :: i, ii

c  Do point product & sum:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = cz*z(i) 
          z(i) = x(i) * y(i) + prod
        enddo
      enddo
#endif
      do i = 1, nxy
        z(i) = x(i) * y(i) + cz*z(i)
      enddo

      return
c
c  end of subroutine qvpvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qdotg
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qdotg(rdot,x,y,nxy,isz) bind(C,name="qdotg")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent (in) :: isz
      integer  (C_SIZE_T), intent (in) :: nxy
      real(C_LONG_DOUBLE), intent (in) :: x(*), y(*)
      real(C_LONG_DOUBLE), intent(out) :: rdot
      real(C_LONG_DOUBLE)              :: e, s, yy, z
      integer  (C_SIZE_T)              :: i, ii

c  Do vector dot product:
      rdot = 0.0
      e    = 0.0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          yy   = x(i)*y(i) - e
          z    = rdot + yy
          e    = (z - rdot) - yy
          rdot = z
        enddo
      enddo

      return
c
c  end of subroutine qdotg
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qcopy
c DESCRIPTION: copies y into x
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qcopy(x, y, nxy, isz) bind(C,name="qcopy")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer     (C_INT), intent (in) :: isz
      integer  (C_SIZE_T), intent (in) :: nxy
      real(C_LONG_DOUBLE), intent (in) :: y(*)
      real(C_LONG_DOUBLE), intent(out) :: x(*)
      integer  (C_SIZE_T)              :: i, ii

c  Do vector copy:
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          x(i) = y(i) 
        enddo
      enddo

      return
c
c  end of subroutine qcopy
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qoop, for real*8
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qoop(u,g,ig,nm,nop) bind(C,name="qoop")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer  (C_SIZE_T), intent (in) :: ig(nop*nm), nm, nop
      real(C_LONG_DOUBLE), intent (in) :: g (*)
      real(C_LONG_DOUBLE), intent(out) :: u (*)
      real(C_LONG_DOUBLE)              :: gtmp
      integer  (C_SIZE_T)              :: i, j, k
c
c     stop "subroutine qoop: do not go there!"
c  Do assignment via indirection:  
c     write(*,*), "nm==", nm
c     write(*,*), "nop=", nop
c     isz =  1
      do i = 1, nop
        gtmp = u(i)
        k    = nm*i
c     write(*,*), "g=", gtmp
c       nn   = im(i)
        do j = 1, nm
c     write(*,*), "ig=", ig(k+j)
          u(ig(k+j)+1) = gtmp
c     write(*,*), "u=", u(ig(i*nm+j))
        enddo
      enddo
c
      return
c
c  end of subroutine qoop
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qisassign
c DESCRIPTION: do assignment by strided-indirection, for real*8
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qisassign(u,g,ig,nm,nd) bind(C,name="qisassign")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer  (C_SIZE_T), intent (in) :: nd, nm, ig(*)
      real(C_LONG_DOUBLE), intent (in) :: g (*)
      real(C_LONG_DOUBLE), intent(out) :: u (*)
      real(C_LONG_DOUBLE)              :: gtmp
      integer  (C_SIZE_T)              :: i, j, k, m
c
c  Do assignment via indirection:  
c     isz =  1
      do i = 1, nd
        gtmp = u(i)
        k    = nm*i
c       nn   = im(i)
        do j = 1, nm  ! stride for this i
          m = ig(k+j)+1
          u(m) = gtmp
        enddo
      enddo
c
      return
c
c  end of subroutine qisassign
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qisum
c DESCRIPTION: do addition using indirection, or real*8
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qisum(res, u, ig, ne) bind(C,name="qisum")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer  (C_SIZE_T), intent (in) :: ne, ig(*)
      real(C_LONG_DOUBLE), intent (in) :: u (*)
      real(C_LONG_DOUBLE), intent(out) :: res(*)
      integer  (C_SIZE_T)              :: j
c
c  Do addition via indirection:  
c     isz =  1
      do j = 1, ne
        res(j) = 0.0;
      enddo
      do j = 1, ne
        res(j) = res(j) + u(ig(j)+1) 
      enddo
c
      return
c
c  end of subroutine qisum
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qiprod
c DESCRIPTION: do multiplication using indirection, for real*8 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qiprod(res, u, ig, ne) bind(C,name="qiprod")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer  (C_SIZE_T), intent (in) :: ne, ig(*)
      real(C_LONG_DOUBLE), intent (in) :: u (*)
      real(C_LONG_DOUBLE), intent(out) :: res
      integer  (C_SIZE_T)              :: j
c
c  Do addition via indirection:  
c     isz =  1
      res = 1.0;
      do j = 1, ne
        res = res * u(ig(j)+1)
      enddo
c
      return
c
c  end of subroutine qiprod
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qimax
c DESCRIPTION: find maxiumum check using indirection 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qimax(res, u, ig, ne) bind(C,name="qimax")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer  (C_SIZE_T), intent (in) :: ne, ig(*)
      real(C_LONG_DOUBLE), intent (in) :: u (*)
      real(C_LONG_DOUBLE), intent(out) :: res
      integer  (C_SIZE_T)              :: j
c
c  Do max via indirection:  
c     isz =  1
      res = -1.0e15
      do j = 1, ne
        res = MAX( res, u(ig(j)+1) )
      enddo
c
      return
c
c  end of subroutine qimax
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : qimin
c DESCRIPTION: find minimum check using indirection, for real*8 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qimin(res, u, ig, ne) bind(C,name="qimin")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer  (C_SIZE_T), intent (in) :: ne, ig(*)
      real(C_LONG_DOUBLE), intent (in) :: u (*)
      real(C_LONG_DOUBLE), intent(out) :: res
      integer  (C_SIZE_T)              :: j
c
c  Do min via indirection:  
c     isz =  1
      res = 1.0e15
      do j = 1, ne
        res = MIN( res, u(ig(j)+1) )
      enddo
c
      return
c
c  end of subroutine qimin
c
      end


c======================================================================================
c Double-precision versions
c======================================================================================


c************************************************************************************
c************************************************************************************
c METHOD     : dmxm
c DESCRIPTION: mxm driver routine that that calls methods explicitly unrolled
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm(C,A,nai,naj,B,nbi,nbj,isz) bind(C,name="dmxm")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      integer   (C_INT), intent (in) :: isz
      integer(C_SIZE_T), intent (in) :: nai, naj, nbi, nbj
      real   (C_DOUBLE), intent (in) :: A(nai,naj), B(nbi,nbj)
      real   (C_DOUBLE), intent(out) :: C(nai,nbj)
      real   (C_DOUBLE)              :: alpha, beta

!!$acc routine vector dmxm1
!!$acc routine vector dmxm2
!!$acc routine vector dmxm3
!!$acc routine vector dmxm4
!!$acc routine vector dmxm5
!!$acc routine vector dmxm6
!!$acc routine vector dmxm7
!!$acc routine vector dmxm8
!!$acc routine vector dmxm9
!!$acc routine vector dmxm10
!!$acc routine vector dmxm11
!!$acc routine vector dmxm12
!!$acc routine vector dmxm13
!!$acc routine vector dmxm14
!!$acc routine vector dmxm15
!!$acc routine vector dmxm16
!!$acc routine vector dmxm17
!!$acc routine vector dmxm18
!!$acc routine vector dmxm19
!!$acc routine vector dmxm20
!!$acc routine vector dmxm21
!!$acc routine vector dmxm22
!!$acc routine vector dmxm23
!!$acc routine vector dmxm24
!!!$acc routine vector dmxmcf
      if      ( nbi .eq.  1 ) then
        call dmxm1 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  2 ) then
        call dmxm2 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  3 ) then
        call dmxm3 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  4 ) then
        call dmxm4 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  5 ) then
        call dmxm5 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  6 ) then
        call dmxm6 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  7 ) then
        call dmxm7 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  8 ) then
        call dmxm8 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  9 ) then
        call dmxm9 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 10 ) then
        call dmxm10(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 11 ) then
        call dmxm11(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 12 ) then
        call dmxm12(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 13 ) then
        call dmxm13(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 14 ) then
        call dmxm14(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 15 ) then
        call dmxm15(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 16 ) then
        call dmxm16(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 17 ) then
        call dmxm17(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 18 ) then
        call dmxm18(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 19 ) then
        call dmxm19(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 20 ) then
        call dmxm20(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 21 ) then
        call dmxm21(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 22 ) then
        call dmxm22(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 23 ) then
        call dmxm23(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 24 ) then
        call dmxm24(A, nai, B, nbj, C, nbj)
      else
        call dmxmcf(C, A, nai, naj, B, nbi, nbj, isz )
#if 0
        alpha = 1.0
        beta  = 0.0
        call dgemm("N","N", nai, nbj, alpha, A, nai, 
     1             B, nbi, beta, C, nai)
#endif
      endif


      return
c
c  end of subroutine dmxm
c
      end subroutine dmxm


c************************************************************************************
c************************************************************************************
c METHOD     : dmxmcf
c DESCRIPTION: Pure Fortran cache-friendly mxm routine that does not unroll by hand,
c              but instead uses cache blocking for mat-mat multiply
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxmcf(C,A,nai,naj,B,nbi,nbj,isz)
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent   (in) :: nai, naj, nbi, nbj
      integer   (C_INT), intent   (in) :: isz
      real   (C_DOUBLE), intent   (in) :: A(nai,naj), B(nbi,nbj)
      real   (C_DOUBLE), intent  (out) :: C(nai,nbj)
      real   (C_DOUBLE)                :: zero
      integer(C_SIZE_T)                :: i, ii, j, jj, k, kk
!     integer*8, intent   (in) :: nai, naj, nbi, nbj, isz
!     real*8   , intent   (in) :: A(nai,naj), B(nbi,nbj)
!     real*8   , intent  (out) :: C(nai,nbj)
!     real*8                   :: zero
!     integer*8                :: i, ii, j, jj, k, kk

c  Initialize C to 0:
      !$acc parallel

!     zero = 0.0
!     call matset(C, nai, nbj, zero, isz)
!     'matset' seems to be required to be provided here 
!     explicitly as ACC isn't finding its device code:
      do ii = 1, nai, isz
        do jj = 1, nbj, isz

          do i = ii, min(nai,ii+isz-1)
            do j = jj, min(nbj,jj+isz-1)
                C(i,j) = 0.0
            enddo
          enddo

        enddo
      enddo


c  Do matrix-matrix multiply using cache-blocking:
      do ii = 1, nai, isz
        do jj = 1, nbj, isz
          do kk = 1, nbi, isz
  
            do i = ii, min(nai,ii+isz-1)
              do j = jj, min(nbj,jj+isz-1)
                do k = kk, min(nbi,kk+isz-1)
                  C(i,j) = C(i,j) + A(i,k)*B(k,j)
                enddo
              enddo
            enddo
        
          enddo
        enddo
      enddo

      !$acc end parallel

      return
c
c  end of subroutine dmxmcf
c
      end subroutine dmxmcf

c************************************************************************************
c************************************************************************************
c METHOD     : dmxv
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxv(y, A, x, nx, ny, isz)  bind(C,name="dmxv")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent  (in)     :: isz
      integer(C_SIZE_T), intent  (in)     :: ny, nx
      real   (C_DOUBLE), intent  (in)     :: A(ny,*), x(*)
      real   (C_DOUBLE), intent (out)     :: y(*)
      integer(C_SIZE_T)                   :: i, ii, j, jj

      do ii = 1, ny, isz
        do i = ii, min(ny,ii+isz-1)
           y(i) = 0.0
        enddo
      enddo


c  Do matrix-vector multiply:
      do ii = 1, ny, isz
        do jj = 1, nx, isz

          do i = ii, min(ny,ii+isz-1)
            do j = jj, min(nx,jj+isz-1)
              y(i) = y(i) + A(i,j)*x(j)
            enddo
          enddo

        enddo
      enddo

      return
c
c  end of subroutine dmxv
c
      end subroutine dmxv

c************************************************************************************
c************************************************************************************
c METHOD     : dmxDm  
c DESCRIPTION: Computes C = A * Diag(B). Dimensions of C assumed to be 
c              those of A
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxDm(C,A,nai,naj,b,nb,isz) bind(C,name="dmxDm")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in)  ::  isz
      integer(C_SIZE_T), intent (in)  ::  nai, naj, nb
      real   (C_DOUBLE), intent (in)  ::  A(nai,naj), b(*)
      real   (C_DOUBLE), intent(out)  ::  C(nai,naj)

      integer(C_SIZE_T)               ::  i, ii, j, jj

c  Do matrix-Diag(matrix) multiply:
      do ii = 1, nai, isz 
        do jj = 1, naj, isz 

          do i = ii, min(nai,ii+isz-1)
            do j = jj, min(naj,jj+isz-1)
              C(i,j) = A(i,j)*b(j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine dmxDm
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : dDmxm  
c DESCRIPTION: Computes C = Diag(A)*B. Dimensions of C assumed to be 
c              those of B
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dDmxm(C, a, na, B, nbi, nbj, isz) bind(C,name="dDmxm")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) ::  isz
      integer(C_SIZE_T), intent (in) ::  nbi, nbj, na
      real   (C_DOUBLE), intent (in) ::  B(nbi,nbj), a(na)
      real   (C_DOUBLE), intent(out) ::  C(nbi,nbj)
      integer(C_SIZE_T)              ::  i, ii, j, jj

c  Do Diag(matrix)-matrix multiply:
      do ii = 1, nbi, isz 
        do jj = 1, nbj, isz 

          do i = ii, min(nbi,ii+isz-1)
            do j = jj, min(nbj,jj+isz-1)
              C(i,j) = a(i)*B(i,j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine dDmxm
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : daApbB  
c DESCRIPTION: computes C = a A + b B. Dimensions of A, B, C assumed to be equal
c              (nxm), and a, and b are double constants
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine daApbB(C,A,B,n,m,ca,cb,isz) bind(C,name="daApbB")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) :: isz
      integer(C_SIZE_T), intent (in) :: n, m
      real   (C_DOUBLE), intent (in) :: A(n,m), B(n,m), ca, cb
      real   (C_DOUBLE), intent(out) :: C(n,m)
      integer(C_SIZE_T)              :: i, ii, j, jj

c  Do sum:
      do ii = 1, m, isz 
        do jj = 1, n, isz 

          do i = ii, min(m,ii+isz-1)
            do j = jj, min(n,jj+isz-1)
              C(i,j) = ca*A(i,j) + cb*B(i,j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine daApbB
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : dzaxpby
c DESCRIPTION: z = a*x + b*y
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dzaxpby(z,x,a,y,b,nxy,isz) bind(C,name="dzaxpby")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) ::  isz
      integer(C_SIZE_T), intent (in) ::  nxy
      real   (C_DOUBLE), intent (in) :: x(nxy), y(nxy), a, b
      real   (C_DOUBLE), intent(out) :: z(nxy)
      real   (C_DOUBLE)              :: prod
      integer(C_SIZE_T)              :: i, ii
c
c  Do vector sums:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = b*y(i) 
          z(i) = a*x(i) + prod
        enddo
      enddo
#endif

      do i = 1, nxy
        z(i) = a*x(i) + b*y(i)
      enddo

      return
c
c  end of subroutine dzaxpby
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dxaxpby
c DESCRIPTION: x(i) = a*x(i)+b*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dxaxpby(x,a,y,b,nxy,isz) bind(C,name="dxaxpby")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent   (in) :: isz
      integer(C_SIZE_T), intent   (in) :: nxy
      real   (C_DOUBLE), intent   (in) :: y(*), a, b
      real   (C_DOUBLE), intent(inout) :: x(*)
      real   (C_DOUBLE)                :: prod
      integer(C_SIZE_T)                :: i, ii

c  Do vector sums:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = b*y(i) + 0.0
          x(i) = a*x(i) + prod
        enddo
      enddo
#endif
      do i = 1, nxy
        x(i) = a*x(i) + b*y(i)
      enddo

      return
c
c  end of subroutine dxaxpby
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dzvxvpt
c DESCRIPTION: z(i) = x(i)*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dzvxvpt(z,x,y,nxy,isz) bind(C,name="dzvxvpt")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) :: isz
      integer(C_SIZE_T), intent (in) :: nxy
      real   (C_DOUBLE), intent (in) :: x(*), y(*)
      real   (C_DOUBLE), intent(out) :: z(*)
      integer(C_SIZE_T)              :: i, ii

c  Do point product:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          z(i) = x(i) * y(i) 
        enddo
      enddo
#endif

      do i = 1, nxy
        z(i) = x(i) * y(i) 
      enddo

      return
c
c  end of subroutine dzvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dvvxvpt
c DESCRIPTION: x(i) = x(i)*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dvvxvpt(x,y,nxy,isz) bind(C,name="dvvxvpt")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent   (in) :: isz
      integer(C_SIZE_T), intent   (in) :: nxy
      real   (C_DOUBLE), intent   (in) :: y(*)
      real   (C_DOUBLE), intent(inout) :: x(*)
      integer(C_SIZE_T)                :: i, ii

c  Do point product:
#if 0
      do ii = 1, nxy, isz

        do i = ii, min(nxy,ii+isz-1)
          x(i) = x(i) * y(i) 
        enddo

      enddo
#endif

      do i = 1, nxy
        x(i) = x(i) * y(i) 
      enddo

      return
c
c  end of subroutine dvvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dvvxvptpv
c DESCRIPTION: x(i) = x(i)*y(i) + cz*z(i)
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dvvxvptpv(x,y,z,cz,nxy,isz) bind(C,name="dvvxvptpv")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent   (in) :: isz
      integer(C_SIZE_T), intent   (in) :: nxy
      real   (C_DOUBLE), intent   (in) :: y(*), z(*), cz
      real   (C_DOUBLE), intent(inout) :: x(*)
      real   (C_DOUBLE)                :: prod
      integer(C_SIZE_T)                :: i, ii

c  Do point product & sum:
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = cz*z(i) 
          x(i) = x(i) * y(i) + prod
        enddo
      enddo

      return
c
c  end of subroutine dvvxvptpv
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dvpvxvpt
c DESCRIPTION: z(i) = x(i)*y(i) + cz*z(i)
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dvpvxvpt(x,y,z,cz,nxy,isz) bind(C,name="dvpvxvpt")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent   (in) :: isz
      integer(C_SIZE_T), intent   (in) :: nxy
      real   (C_DOUBLE), intent   (in) :: y(*), x(*), cz
      real   (C_DOUBLE), intent(inout) :: z(*)
      real   (C_DOUBLE)                :: prod
      integer(C_SIZE_T)                :: i, ii

c  Do point product & sum:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = cz*z(i) 
          z(i) = x(i) * y(i) + prod
        enddo
      enddo
#endif
      do i = 1, nxy
        z(i) = x(i) * y(i) + cz*z(i)
      enddo

      return
c
c  end of subroutine dvpvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : ddotg
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine ddotg(rdot,x,y,nxy,isz) bind(C,name="ddotg")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) :: isz
      integer(C_SIZE_T), intent (in) :: nxy
      real   (C_DOUBLE), intent (in) :: x(*), y(*)
      real   (C_DOUBLE), intent(out) :: rdot
      real   (C_DOUBLE)              :: e, s, yy, z
      integer(C_SIZE_T)              :: i, ii

c  Do vector dot product:
      rdot = 0.0
      e    = 0.0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          yy   = x(i)*y(i) - e
          z    = rdot + yy
          e    = (z - rdot) - yy
          rdot = z
        enddo
      enddo

      return
c
c  end of subroutine ddotg
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dcopy
c DESCRIPTION: copies y into x
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dcopy(x, y, nxy, isz) bind(C,name="dcopy")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) :: isz
      integer(C_SIZE_T), intent (in) :: nxy
      real   (C_DOUBLE), intent (in) :: y(*)
      real   (C_DOUBLE), intent(out) :: x(*)
      integer(C_SIZE_T)              :: i, ii

c  Do vector copy:
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          x(i) = y(i) 
        enddo
      enddo

      return
c
c  end of subroutine dcopy
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : matset
c DESCRIPTION: sets matrix to specified value
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine matset(C,nci,ncj,val,isz) 
      implicit none
!$acc routine vector
      integer , intent (in) :: nci, ncj, isz
      real    , intent (in) :: val
      real    , intent(out) :: C(nci,ncj)
      integer               :: i, j, ii, jj

      do ii = 1, nci, isz
        do jj = 1, ncj, isz

          do i = ii, min(nci,ii+isz-1)
            do j = jj, min(ncj,jj+isz-1)
                C(i,j) = val
            enddo
          enddo

        enddo
      enddo

      return
c
c  end of subroutine matset
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : doop, for real*8
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine doop(u,g,ig,nm,nop) bind(C,name="doop")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: ig(nop*nm), nm, nop
      real   (C_DOUBLE), intent (in) :: g (*)
      real   (C_DOUBLE), intent(out) :: u (*)
      real   (C_DOUBLE)              :: gtmp
      integer(C_SIZE_T)              :: i, j, k
c
c     stop "subroutine doop: do not go there!"
c  Do assignment via indirection:  
c     write(*,*), "nm==", nm
c     write(*,*), "nop=", nop
c     isz =  1
      do i = 1, nop
        gtmp = u(i)
        k    = nm*i
c     write(*,*), "g=", gtmp
c       nn   = im(i)
        do j = 1, nm
c     write(*,*), "ig=", ig(k+j)
          u(ig(k+j)+1) = gtmp
c     write(*,*), "u=", u(ig(i*nm+j))
        enddo
      enddo
c
      return
c
c  end of subroutine doop
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : disassign
c DESCRIPTION: do assignment by strided-indirection, for real*8
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine disassign(u,g,ig,nm,nd) bind(C,name="disassign")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: nd, nm, ig(*)
      real   (C_DOUBLE), intent (in) :: g (*)
      real   (C_DOUBLE), intent(out) :: u (*)
      real   (C_DOUBLE)              :: gtmp
      integer(C_SIZE_T)              :: i, j, k, m
c
c  Do assignment via indirection:  
c     isz =  1
      do i = 1, nd
        gtmp = u(i)
        k    = nm*i
c       nn   = im(i)
        do j = 1, nm  ! stride for this i
          m = ig(k+j)+1
          u(m) = gtmp
        enddo
      enddo
c
      return
c
c  end of subroutine disassign
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : disum
c DESCRIPTION: do addition using indirection, or real*8
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine disum(res, u, ig, ne) bind(C,name="disum")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: ne, ig(*)
      real   (C_DOUBLE), intent (in) :: u (*)
      real   (C_DOUBLE), intent(out) :: res(*)
      integer(C_SIZE_T)              :: j
c
c  Do addition via indirection:  
c     isz =  1
      do j = 1, ne
        res(j) = 0.0;
      enddo
      do j = 1, ne
        res(j) = res(j) + u(ig(j)+1) 
      enddo
c
      return
c
c  end of subroutine disum
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : diprod
c DESCRIPTION: do multiplication using indirection, for real*8 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine diprod(res, u, ig, ne) bind(C,name="diprod")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: ne, ig(*)
      real   (C_DOUBLE), intent (in) :: u (*)
      real   (C_DOUBLE), intent(out) :: res
      integer(C_SIZE_T)              :: j
c
c  Do addition via indirection:  
c     isz =  1
      res = 1.0;
      do j = 1, ne
        res = res * u(ig(j)+1)
      enddo
c
      return
c
c  end of subroutine diprod
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : imax
c DESCRIPTION: find maxiumum check using indirection 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dimax(res, u, ig, ne) bind(C,name="imax")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: ne, ig(*)
      real   (C_DOUBLE), intent (in) :: u (*)
      real   (C_DOUBLE), intent(out) :: res
      integer(C_SIZE_T)              :: j
c
c  Do max via indirection:  
c     isz =  1
      res = -1.0e15
      do j = 1, ne
        res = MAX( res, u(ig(j)+1) )
      enddo
c
      return
c
c  end of subroutine imax
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dimin
c DESCRIPTION: find minimum check using indirection, for real*8 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dimin(res, u, ig, ne) bind(C,name="dimin")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: ne, ig(*)
      real   (C_DOUBLE), intent (in) :: u (*)
      real   (C_DOUBLE), intent(out) :: res
      integer(C_SIZE_T)              :: j
c
c  Do min via indirection:  
c     isz =  1
      res = 1.0e15
      do j = 1, ne
        res = MIN( res, u(ig(j)+1) )
      enddo
c
      return
c
c  end of subroutine dimin
c
      end


c======================================================================================
c Float/ single precision versions
c======================================================================================



c************************************************************************************
c************************************************************************************
c METHOD     : fmxm
c DESCRIPTION: mxm driver routine that that calls methods explicitly unrolled
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm(C,A,nai,naj,B,nbi,nbj,isz) bind(C,name="fmxm")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      integer   (C_INT), intent (in) :: isz
      integer(C_SIZE_T), intent (in) :: nai, naj, nbi, nbj
      real    (C_FLOAT), intent (in) :: A(nai,naj), B(nbi,nbj)
      real    (C_FLOAT), intent(out) :: C(nai,nbj)
      real    (C_FLOAT)              :: alpha, beta

!!$acc routine vector fmxm1
!!$acc routine vector fmxm2
!!$acc routine vector fmxm3
!!$acc routine vector fmxm4
!!$acc routine vector fmxm5
!!$acc routine vector fmxm6
!!$acc routine vector fmxm7
!!$acc routine vector fmxm8
!!$acc routine vector fmxm9
!!$acc routine vector fmxm10
!!$acc routine vector fmxm11
!!$acc routine vector fmxm12
!!$acc routine vector fmxm13
!!$acc routine vector fmxm14
!!$acc routine vector fmxm15
!!$acc routine vector fmxm16
!!$acc routine vector fmxm17
!!$acc routine vector fmxm18
!!$acc routine vector fmxm19
!!$acc routine vector fmxm20
!!$acc routine vector fmxm21
!!$acc routine vector fmxm22
!!$acc routine vector fmxm23
!!$acc routine vector fmxm24
!!!$acc routine vector fmxmcf
      if      ( nbi .eq.  1 ) then
        call fmxm1 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  2 ) then
        call fmxm2 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  3 ) then
        call fmxm3 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  4 ) then
        call fmxm4 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  5 ) then
        call fmxm5 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  6 ) then
        call fmxm6 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  7 ) then
        call fmxm7 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  8 ) then
        call fmxm8 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  9 ) then
        call fmxm9 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 10 ) then
        call fmxm10(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 11 ) then
        call fmxm11(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 12 ) then
        call fmxm12(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 13 ) then
        call fmxm13(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 14 ) then
        call fmxm14(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 15 ) then
        call fmxm15(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 16 ) then
        call fmxm16(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 17 ) then
        call fmxm17(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 18 ) then
        call fmxm18(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 19 ) then
        call fmxm19(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 20 ) then
        call fmxm20(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 21 ) then
        call fmxm21(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 22 ) then
        call fmxm22(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 23 ) then
        call fmxm23(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 24 ) then
        call fmxm24(A, nai, B, nbj, C, nbj)
      else
        call fmxmcf(C, A, nai, naj, B, nbi, nbj, isz )
#if 0
        alpha = 1.0
        beta  = 0.0
        call dgemm("N","N", nai, nbj, alpha, A, nai, 
     1             B, nbi, beta, C, nai)
#endif
      endif


      return
c
c  end of subroutine fmxm
c
      end subroutine fmxm


c************************************************************************************
c************************************************************************************
c METHOD     : fmxmcf
c DESCRIPTION: Pure Fortran cache-friendly mxm routine that does not unroll by hand,
c              but instead uses cache blocking for mat-mat multiply
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxmcf(C,A,nai,naj,B,nbi,nbj,isz) 
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent   (in) :: nai, naj, nbi, nbj
      integer   (C_INT), intent   (in) :: isz
      real    (C_FLOAT), intent   (in) :: A(nai,naj), B(nbi,nbj)
      real    (C_FLOAT), intent  (out) :: C(nai,nbj)
      real    (C_FLOAT)                :: zero
      integer(C_SIZE_T)                :: i, ii, j, jj, k, kk
!     integer*8, intent   (in) :: nai, naj, nbi, nbj, isz
!     real*8   , intent   (in) :: A(nai,naj), B(nbi,nbj)
!     real*8   , intent  (out) :: C(nai,nbj)
!     real*8                   :: zero
!     integer*8                :: i, ii, j, jj, k, kk

c  Initialize C to 0:
      !$acc parallel

!     zero = 0.0
!     call matset(C, nai, nbj, zero, isz)
!     'matset' seems to be required to be provided here 
!     explicitly as ACC isn't finding its device code:
      do ii = 1, nai, isz
        do jj = 1, nbj, isz

          do i = ii, min(nai,ii+isz-1)
            do j = jj, min(nbj,jj+isz-1)
                C(i,j) = 0.0
            enddo
          enddo

        enddo
      enddo


c  Do matrix-matrix multiply using cache-blocking:
      do ii = 1, nai, isz
        do jj = 1, nbj, isz
          do kk = 1, nbi, isz
  
            do i = ii, min(nai,ii+isz-1)
              do j = jj, min(nbj,jj+isz-1)
                do k = kk, min(nbi,kk+isz-1)
                  C(i,j) = C(i,j) + A(i,k)*B(k,j)
                enddo
              enddo
            enddo
        
          enddo
        enddo
      enddo

      !$acc end parallel

      return
c
c  end of subroutine fmxmcf
c
      end subroutine fmxmcf

c************************************************************************************
c************************************************************************************
c METHOD     : fmxv
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxv(y, A, x, nx, ny, isz)  bind(C,name="fmxv")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent  (in)     :: isz
      integer(C_SIZE_T), intent  (in)     :: ny, nx
      real    (C_FLOAT), intent  (in)     :: A(ny,*), x(*)
      real    (C_FLOAT), intent (out)     :: y(*)
      integer(C_SIZE_T)                   :: i, ii, j, jj

      do ii = 1, ny, isz
        do i = ii, min(ny,ii+isz-1)
           y(i) = 0.0
        enddo
      enddo


c  Do matrix-vector multiply:
      do ii = 1, ny, isz
        do jj = 1, nx, isz

          do i = ii, min(ny,ii+isz-1)
            do j = jj, min(nx,jj+isz-1)
              y(i) = y(i) + A(i,j)*x(j)
            enddo
          enddo

        enddo
      enddo

      return
c
c  end of subroutine fmxv
c
      end subroutine fmxv

c************************************************************************************
c************************************************************************************
c METHOD     : fmxDm  
c DESCRIPTION: computes C = A * Diag(B). Dimensions of C assumed to be 
c              those of A
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxDm(C,A,nai,naj,b,nb,isz) bind(C,name="fmxDm")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in)  ::  isz
      integer(C_SIZE_T), intent (in)  ::  nai, naj, nb
      real    (C_FLOAT), intent (in)  ::  A(nai,naj), b(*)
      real    (C_FLOAT), intent(out)  ::  C(nai,naj)

      integer(C_SIZE_T)               ::  i, ii, j, jj

c  Do matrix-Diag(matrix) multiply:
      do ii = 1, nai, isz 
        do jj = 1, naj, isz 

          do i = ii, min(nai,ii+isz-1)
            do j = jj, min(naj,jj+isz-1)
              C(i,j) = A(i,j)*b(j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine fmxDm
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : fDmxm  
c DESCRIPTION: computes C = Diag(A)*B. Dimensions of C assumed to be 
c              those of B
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fDmxm(C, a, na, B, nbi, nbj, isz) bind(C,name="fDmxm")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) ::  isz
      integer(C_SIZE_T), intent (in) ::  nbi, nbj, na
      real    (C_FLOAT), intent (in) ::  B(nbi,nbj), a(na)
      real    (C_FLOAT), intent(out) ::  C(nbi,nbj)
      integer(C_SIZE_T)              ::  i, ii, j, jj

c  Do Diag(matrix)-matrix multiply:
      do ii = 1, nbi, isz 
        do jj = 1, nbj, isz 

          do i = ii, min(nbi,ii+isz-1)
            do j = jj, min(nbj,jj+isz-1)
              C(i,j) = a(i)*B(i,j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine fDmxm
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : faApbB  
c DESCRIPTION: computes C = a A + b B. Dimensions of A, B, C assumed to be equal
c              (nxm), and a, and b are double constants
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine faApbB(C,A,B,n,m,ca,cb,isz) bind(C,name="faApbB")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) :: isz
      integer(C_SIZE_T), intent (in) :: n, m
      real    (C_FLOAT), intent (in) :: A(n,m), B(n,m), ca, cb
      real    (C_FLOAT), intent(out) :: C(n,m)
      integer(C_SIZE_T)              :: i, ii, j, jj

c  Do sum:
      do ii = 1, m, isz 
        do jj = 1, n, isz 

          do i = ii, min(m,ii+isz-1)
            do j = jj, min(n,jj+isz-1)
              C(i,j) = ca*A(i,j) + cb*B(i,j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine faApbB
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : fzaxpby
c DESCRIPTION: z = a*x + b*y
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fzaxpby(z,x,a,y,b,nxy,isz) bind(C,name="fzaxpby")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) ::  isz
      integer(C_SIZE_T), intent (in) ::  nxy
      real    (C_FLOAT), intent (in) :: x(nxy), y(nxy), a, b
      real    (C_FLOAT), intent(out) :: z(nxy)
      real    (C_FLOAT)              :: prod
      integer(C_SIZE_T)              :: i, ii
c
c  Do vector sums:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = b*y(i) 
          z(i) = a*x(i) + prod
        enddo
      enddo
#endif

      do i = 1, nxy
        z(i) = a*x(i) + b*y(i)
      enddo

      return
c
c  end of subroutine fzaxpby
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fxaxpby
c DESCRIPTION: x(i) = a*x(i)+b*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fxaxpby(x,a,y,b,nxy,isz) bind(C,name="fxaxpby")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent   (in) :: isz
      integer(C_SIZE_T), intent   (in) :: nxy
      real    (C_FLOAT), intent   (in) :: y(*), a, b
      real    (C_FLOAT), intent(inout) :: x(*)
      real    (C_FLOAT)                :: prod
      integer(C_SIZE_T)                :: i, ii

c  Do vector sums:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = b*y(i) + 0.0
          x(i) = a*x(i) + prod
        enddo
      enddo
#endif
      do i = 1, nxy
        x(i) = a*x(i) + b*y(i)
      enddo

      return
c
c  end of subroutine fxaxpby
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fzvxvpt
c DESCRIPTION: z(i) = x(i)*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fzvxvpt(z,x,y,nxy,isz) bind(C,name="fzvxvpt")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) :: isz
      integer(C_SIZE_T), intent (in) :: nxy
      real    (C_FLOAT), intent (in) :: x(*), y(*)
      real    (C_FLOAT), intent(out) :: z(*)
      integer(C_SIZE_T)              :: i, ii

c  Do point product:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          z(i) = x(i) * y(i) 
        enddo
      enddo
#endif

      do i = 1, nxy
        z(i) = x(i) * y(i) 
      enddo

      return
c
c  end of subroutine fzvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fvvxvpt
c DESCRIPTION: x(i) = x(i)*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fvvxvpt(x,y,nxy,isz) bind(C,name="fvvxvpt")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent   (in) :: isz
      integer(C_SIZE_T), intent   (in) :: nxy
      real    (C_FLOAT), intent   (in) :: y(*)
      real    (C_FLOAT), intent(inout) :: x(*)
      integer(C_SIZE_T)                :: i, ii

c  Do point product:
#if 0
      do ii = 1, nxy, isz

        do i = ii, min(nxy,ii+isz-1)
          x(i) = x(i) * y(i) 
        enddo

      enddo
#endif

      do i = 1, nxy
        x(i) = x(i) * y(i) 
      enddo

      return
c
c  end of subroutine fvvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fvvxvptpv
c DESCRIPTION: x(i) = x(i)*y(i) + cz*z(i)
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fvvxvptpv(x,y,z,cz,nxy,isz) bind(C,name="fvvxvptpv")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent   (in) :: isz
      integer(C_SIZE_T), intent   (in) :: nxy
      real    (C_FLOAT), intent   (in) :: y(*), z(*), cz
      real    (C_FLOAT), intent(inout) :: x(*)
      real    (C_FLOAT)                :: prod
      integer(C_SIZE_T)                :: i, ii

c  Do point product & sum:
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = cz*z(i) 
          x(i) = x(i) * y(i) + prod
        enddo
      enddo

      return
c
c  end of subroutine fvvxvptpv
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fvpvxvpt
c DESCRIPTION: z(i) = x(i)*y(i) + cz*z(i)
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fvpvxvpt(x,y,z,cz,nxy,isz) bind(C,name="fvpvxvpt")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent   (in) :: isz
      integer(C_SIZE_T), intent   (in) :: nxy
      real    (C_FLOAT), intent   (in) :: y(*), x(*), cz
      real    (C_FLOAT), intent(inout) :: z(*)
      real    (C_FLOAT)                :: prod
      integer(C_SIZE_T)                :: i, ii

c  Do point product & sum:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = cz*z(i) 
          z(i) = x(i) * y(i) + prod
        enddo
      enddo
#endif
      do i = 1, nxy
        z(i) = x(i) * y(i) + cz*z(i)
      enddo

      return
c
c  end of subroutine fvpvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fdotg
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fdotg(rdot,x,y,nxy,isz) bind(C,name="fdotg")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) :: isz
      integer(C_SIZE_T), intent (in) :: nxy
      real    (C_FLOAT), intent (in) :: x(*), y(*)
      real    (C_FLOAT), intent(out) :: rdot
      real    (C_FLOAT)              :: e, s, yy, z
      integer(C_SIZE_T)              :: i, ii

c  Do vector dot product:
      rdot = 0.0
      e    = 0.0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          yy   = x(i)*y(i) - e
          z    = rdot + yy
          e    = (z - rdot) - yy
          rdot = z
        enddo
      enddo

      return
c
c  end of subroutine fdotg
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fcopy
c DESCRIPTION: copies y into x
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fcopy(x, y, nxy, isz) bind(C,name="fcopy")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer   (C_INT), intent (in) :: isz
      integer(C_SIZE_T), intent (in) :: nxy
      real    (C_FLOAT), intent (in) :: y(*)
      real    (C_FLOAT), intent(out) :: x(*)
      integer(C_SIZE_T)              :: i, ii

c  Do vector copy:
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          x(i) = y(i) 
        enddo
      enddo

      return
c
c  end of subroutine fcopy
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : foop, for real*8
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine foop(u,g,ig,nm,nop) bind(C,name="foop")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: ig(nop*nm), nm, nop
      real    (C_FLOAT), intent (in) :: g (*)
      real    (C_FLOAT), intent(out) :: u (*)
      real    (C_FLOAT)              :: gtmp
      integer(C_SIZE_T)              :: i, j, k
c
c     stop "subroutine foop: do not go there!"
c  Do assignment via indirection:  
c     write(*,*), "nm==", nm
c     write(*,*), "nop=", nop
c     isz =  1
      do i = 1, nop
        gtmp = u(i)
        k    = nm*i
c     write(*,*), "g=", gtmp
c       nn   = im(i)
        do j = 1, nm
c     write(*,*), "ig=", ig(k+j)
          u(ig(k+j)+1) = gtmp
c     write(*,*), "u=", u(ig(i*nm+j))
        enddo
      enddo
c
      return
c
c  end of subroutine foop
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fisassign
c DESCRIPTION: do assignment by strided-indirection, for real*8
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fisassign(u,g,ig,nm,nd) bind(C,name="fisassign")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: nd, nm, ig(*)
      real    (C_FLOAT), intent (in) :: g (*)
      real    (C_FLOAT), intent(out) :: u (*)
      real    (C_FLOAT)              :: gtmp
      integer(C_SIZE_T)              :: i, j, k, m
c
c  Do assignment via indirection:  
c     isz =  1
      do i = 1, nd
        gtmp = u(i)
        k    = nm*i
c       nn   = im(i)
        do j = 1, nm  ! stride for this i
          m = ig(k+j)+1
          u(m) = gtmp
        enddo
      enddo
c
      return
c
c  end of subroutine fisassign
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fisum
c DESCRIPTION: do addition using indirection, or real*8
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fisum(res, u, ig, ne) bind(C,name="fisum")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: ne, ig(*)
      real    (C_FLOAT), intent (in) :: u (*)
      real    (C_FLOAT), intent(out) :: res(*)
      integer(C_SIZE_T)              :: j
c
c  Do addition via indirection:  
c     isz =  1
      do j = 1, ne
        res(j) = 0.0;
      enddo
      do j = 1, ne
        res(j) = res(j) + u(ig(j)+1) 
      enddo
c
      return
c
c  end of subroutine fisum
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fiprod
c DESCRIPTION: do multiplication using indirection, for real*8 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fiprod(res, u, ig, ne) bind(C,name="fiprod")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: ne, ig(*)
      real    (C_FLOAT), intent (in) :: u (*)
      real    (C_FLOAT), intent(out) :: res
      integer(C_SIZE_T)              :: j
c
c  Do addition via indirection:  
c     isz =  1
      res = 1.0;
      do j = 1, ne
        res = res * u(ig(j)+1)
      enddo
c
      return
c
c  end of subroutine fiprod
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fimax
c DESCRIPTION: find maxiumum check using indirection 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fimax(res, u, ig, ne) bind(C,name="fimax")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: ne, ig(*)
      real    (C_FLOAT), intent (in) :: u (*)
      real    (C_FLOAT), intent(out) :: res
      integer(C_SIZE_T)              :: j
c
c  Do max via indirection:  
c     isz =  1
      res = -1.0e15
      do j = 1, ne
        res = MAX( res, u(ig(j)+1) )
      enddo
c
      return
c
c  end of subroutine fimax
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : fimin
c DESCRIPTION: find minimum check using indirection, for real*8 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fimin(res, u, ig, ne) bind(C,name="fimin")
      use, intrinsic :: ISO_C_BINDING
      implicit none
!$acc routine vector
      integer(C_SIZE_T), intent (in) :: ne, ig(*)
      real    (C_FLOAT), intent (in) :: u (*)
      real    (C_FLOAT), intent(out) :: res
      integer(C_SIZE_T)              :: j
c
c  Do min via indirection:  
c     isz =  1
      res = 1.0e15
      do j = 1, ne
        res = MIN( res, u(ig(j)+1) )
      enddo
c
      return
c
c  end of subroutine fimin
c
      end


