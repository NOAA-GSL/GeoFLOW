/*
 * Filename:	spatial.hpp
 * Author:      bflynt
 * Created:		Jan 6, 2021
 * Copyright:   Copyright 2018. Colorado State University. All rights reserved
 */
#ifndef GEOFLOW_TBOX_SPATIAL_H_
#define GEOFLOW_TBOX_SPATIAL_H_

#include "tbox/spatial/bound/box.hpp"
#include "tbox/spatial/bound/concept.hpp"
#include "tbox/spatial/bound/sphere.hpp"

#include "tbox/spatial/common/truncated_multiset.hpp"

#include "tbox/spatial/shared/index/exhaustive.hpp"
#include "tbox/spatial/shared/index/rtree.hpp"

#include "tbox/spatial/shared/index/rtree/algorithm.hpp"
#include "tbox/spatial/shared/index/rtree/leaf.hpp"
#include "tbox/spatial/shared/index/rtree/linear.hpp"
#include "tbox/spatial/shared/index/rtree/node.hpp"
#include "tbox/spatial/shared/index/rtree/page.hpp"
#include "tbox/spatial/shared/index/rtree/quadratic.hpp"

#include "tbox/spatial/shared/predicate/dispatch.hpp"
#include "tbox/spatial/shared/predicate/distance.hpp"
#include "tbox/spatial/shared/predicate/factories.hpp"
#include "tbox/spatial/shared/predicate/spatial.hpp"
#include "tbox/spatial/shared/predicate/tags.hpp"

#endif /* GEOFLOW_TBOX_SPATIAL_H_ */
