/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011-2014, Willow Garage, Inc.
 *  Copyright (c) 2014-2016, Open Source Robotics Foundation
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Open Source Robotics Foundation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/** @author Jia Pan */

#ifndef FCL_TRAVERSAL_MESHSHAPEDISTANCETRAVERSALNODE_H
#define FCL_TRAVERSAL_MESHSHAPEDISTANCETRAVERSALNODE_H

#include "fcl/geometry/shape/compute_bv.h"
#include "fcl/narrowphase/detail/traversal/distance/bvh_shape_distance_traversal_node.h"

namespace fcl
{

namespace detail
{

/// @brief Traversal node for distance between mesh and shape
template <typename BV, typename Shape, typename NarrowPhaseSolver>
class MeshShapeDistanceTraversalNode
    : public BVHShapeDistanceTraversalNode<BV, Shape>
{ 
public:

  using S = typename BV::S;

  MeshShapeDistanceTraversalNode();

  /// @brief Distance testing between leaves (one triangle and one shape)
  void leafTesting(int b1, int b2) const;

  /// @brief Whether the traversal process can stop early
  bool canStop(S c) const;

  Vector3<S>* vertices;
  Triangle* tri_indices;

  S rel_err;
  S abs_err;
    
  const NarrowPhaseSolver* nsolver;
};

template <typename BV, typename Shape, typename NarrowPhaseSolver>
void meshShapeDistanceOrientedNodeLeafTesting(
    int b1, int /* b2 */,
    const BVHModel<BV>* model1,
    const Shape& model2,
    Vector3<typename BV::S>* vertices,
    Triangle* tri_indices,
    const Transform3<typename BV::S>& tf1,
    const Transform3<typename BV::S>& tf2,
    const NarrowPhaseSolver* nsolver,
    bool enable_statistics,
    int & num_leaf_tests,
    const DistanceRequest<typename BV::S>& /* request */,
    DistanceResult<typename BV::S>& result);


template <typename BV, typename Shape, typename NarrowPhaseSolver>
void distancePreprocessOrientedNode(
    const BVHModel<BV>* model1,
    Vector3<typename BV::S>* vertices,
    Triangle* tri_indices,
    int init_tri_id,
    const Shape& model2,
    const Transform3<typename BV::S>& tf1,
    const Transform3<typename BV::S>& tf2,
    const NarrowPhaseSolver* nsolver,
    const DistanceRequest<typename BV::S>& /* request */,
    DistanceResult<typename BV::S>& result);

/// @brief Initialize traversal node for distance computation between one mesh
/// and one shape, given the current transforms
template <typename BV, typename Shape, typename NarrowPhaseSolver>
bool initialize(
    MeshShapeDistanceTraversalNode<BV, Shape, NarrowPhaseSolver>& node,
    BVHModel<BV>& model1,
    Transform3<typename BV::S>& tf1,
    const Shape& model2,
    const Transform3<typename BV::S>& tf2,
    const NarrowPhaseSolver* nsolver,
    const DistanceRequest<typename BV::S>& request,
    DistanceResult<typename BV::S>& result,
    bool use_refit = false, bool refit_bottomup = false);

/// @brief Traversal node for distance between mesh and shape, when mesh BVH is one of the oriented node (RSS, OBBRSS, kIOS)
template <typename Shape, typename NarrowPhaseSolver>
class MeshShapeDistanceTraversalNodeRSS
    : public MeshShapeDistanceTraversalNode<
    RSS<typename Shape::S>, Shape, NarrowPhaseSolver>
{
public:
  using S = typename Shape::S;

  MeshShapeDistanceTraversalNodeRSS();

  void preprocess();

  void postprocess();

  S BVTesting(int b1, int b2) const;

  void leafTesting(int b1, int b2) const;
};

template <typename Shape, typename NarrowPhaseSolver>
bool initialize(
    MeshShapeDistanceTraversalNodeRSS<Shape, NarrowPhaseSolver>& node,
    const BVHModel<RSS<typename Shape::S>>& model1,
    const Transform3<typename Shape::S>& tf1,
    const Shape& model2,
    const Transform3<typename Shape::S>& tf2,
    const NarrowPhaseSolver* nsolver,
    const DistanceRequest<typename Shape::S>& request,
    DistanceResult<typename Shape::S>& result);

template <typename Shape, typename NarrowPhaseSolver>
class MeshShapeDistanceTraversalNodekIOS
    : public MeshShapeDistanceTraversalNode<kIOS<typename Shape::S>, Shape, NarrowPhaseSolver>
{
public:
  using S = typename Shape::S;

  MeshShapeDistanceTraversalNodekIOS();

  void preprocess();

  void postprocess();

  S BVTesting(int b1, int b2) const;

  void leafTesting(int b1, int b2) const;

};

/// @brief Initialize traversal node for distance computation between one mesh and one shape, specialized for kIOS type
template <typename Shape, typename NarrowPhaseSolver>
bool initialize(
    MeshShapeDistanceTraversalNodekIOS<Shape, NarrowPhaseSolver>& node,
    const BVHModel<kIOS<typename Shape::S>>& model1, const Transform3<typename Shape::S>& tf1,
    const Shape& model2, const Transform3<typename Shape::S>& tf2,
    const NarrowPhaseSolver* nsolver,
    const DistanceRequest<typename Shape::S>& request,
    DistanceResult<typename Shape::S>& result);

template <typename Shape, typename NarrowPhaseSolver>
class MeshShapeDistanceTraversalNodeOBBRSS
    : public MeshShapeDistanceTraversalNode<OBBRSS<typename Shape::S>, Shape, NarrowPhaseSolver>
{
public:
  using S = typename Shape::S;

  MeshShapeDistanceTraversalNodeOBBRSS();

  void preprocess();

  void postprocess();

  S BVTesting(int b1, int b2) const;

  void leafTesting(int b1, int b2) const;
  
};

/// @brief Initialize traversal node for distance computation between one mesh and one shape, specialized for OBBRSS type
template <typename Shape, typename NarrowPhaseSolver>
bool initialize(
    MeshShapeDistanceTraversalNodeOBBRSS<Shape, NarrowPhaseSolver>& node,
    const BVHModel<OBBRSS<typename Shape::S>>& model1,
    const Transform3<typename Shape::S>& tf1,
    const Shape& model2,
    const Transform3<typename Shape::S>& tf2,
    const NarrowPhaseSolver* nsolver,
    const DistanceRequest<typename Shape::S>& request,
    DistanceResult<typename Shape::S>& result);

} // namespace detail
} // namespace fcl

#include "fcl/narrowphase/detail/traversal/distance/mesh_shape_distance_traversal_node-inl.h"

#endif
