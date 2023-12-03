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

/** @author Shi Shenglei */

#ifndef FCL_NARROWPHASE_BISECTIONSOLVER_INL_H
#define FCL_NARROWPHASE_BISECTIONSOLVER_INL_H

#include "fcl/narrowphase/detail/bisection_solver.h"

#include <algorithm>

#include "fcl/common/unused.h"

#include "fcl/narrowphase/detail/convexity_based_algorithm/bisection_distance.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/capsule_capsule.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_box.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_capsule.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_cylinder.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_sphere.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_triangle.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/box_box.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/halfspace.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/plane.h"
#include "fcl/narrowphase/detail/failed_at_this_configuration.h"

namespace fcl
{

namespace detail
{

//==============================================================================
extern template
struct BisectionSolver<double>;

//==============================================================================
template<typename S, typename Shape1, typename Shape2>
struct ShapeDistanceBisectionImpl
{
  static bool run(
      const BisectionSolver<S>& gjkSolver,
      const Shape1& s1,
      const Transform3<S>& tf1,
      const Shape2& s2,
      const Transform3<S>& tf2,
      S* dist,
      Vector3<S>* p1,
      Vector3<S>* p2)
  {
    detail::MinkowskiDiff<S> shape;
    shape.shapes[0] = &s1;
    shape.shapes[1] = &s2;
    shape.toshape1.noalias() = tf2.linear().transpose() * tf1.linear();
    shape.toshape0 = tf1.inverse(Eigen::Isometry) * tf2;

    bool res =  detail::SeparationDistance(
          shape,
          gjkSolver.max_distance_iterations,
          gjkSolver.distance_tolerance,
          dist,
          p1,
          p2);

    return res;
  }
};

template<typename S>
template<typename Shape1, typename Shape2>
bool BisectionSolver<S>::shapeDistance(
    const Shape1& s1,
    const Transform3<S>& tf1,
    const Shape2& s2,
    const Transform3<S>& tf2,
    S* dist,
    Vector3<S>* p1,
    Vector3<S>* p2) const
{
  return ShapeDistanceBisectionImpl<S, Shape1, Shape2>::run(
        *this, s1, tf1, s2, tf2, dist, p1, p2);
}

template<typename S>
template<typename Shape1, typename Shape2>
bool BisectionSolver<S>::shapeDistanceCCD(
    const Shape1& s1,
    const Transform3<S>& tf1,
    const Shape2& s2,
    const Transform3<S>& tf2,
    S* dist,
    Vector3<S>* p1,
    Vector3<S>* p2) const
{
  detail::MinkowskiDiff<S> shape;
  shape.shapes[0] = &s1;
  shape.shapes[1] = &s2;
  shape.toshape1.noalias() = tf2.linear().transpose() * tf1.linear();
  shape.toshape0 = tf1.inverse(Eigen::Isometry) * tf2;

  Vector3<S> p1_ = Vector3<S>::Zero(), p2_ = Vector3<S>::Zero();
  S d = detail::libccd_extension::ccdGJKDist2(shape, max_distance_iterations, distance_tolerance, p1_, p2_);
  if (p1) *p1 = p1_;
  if (p2) *p2 = p2_;
  if (dist) *dist = d;
  if (d < 0.0)
    return false;
  else
    return true;
}

// Shape distance algorithms not using libccd
//
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// |            | box | sphere | ellipsoid | capsule | cone | cylinder | plane | half-space | triangle |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | box        |     |   O    |           |         |      |          |       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | sphere     |/////|   O    |           |    O    |      |    O     |       |            |     O    |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | ellipsoid  |/////|////////|           |         |      |          |       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | capsule    |/////|////////|///////////|    O    |      |          |       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | cone       |/////|////////|///////////|/////////|      |          |       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | cylinder   |/////|////////|///////////|/////////|//////|          |       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | plane      |/////|////////|///////////|/////////|//////|//////////|       |            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | half-space |/////|////////|///////////|/////////|//////|//////////|///////|            |          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+
// | triangle   |/////|////////|///////////|/////////|//////|//////////|///////|////////////|          |
// +------------+-----+--------+-----------+---------+------+----------+-------+------------+----------+

//==============================================================================
//template<typename S>
//struct ShapeDistanceBisectionImpl<S, Sphere<S>, Box<S>>
//{
//  static bool run(
//      const BisectionSolver<S>& /*gjkSolver*/,
//      const Sphere<S>& s1,
//      const Transform3<S>& tf1,
//      const Box<S>& s2,
//      const Transform3<S>& tf2,
//      S* dist,
//      Vector3<S>* p1,
//      Vector3<S>* p2)
//  {
//    return detail::sphereBoxDistance(s1, tf1, s2, tf2, dist, p1, p2);
//  }
//};
//
////==============================================================================
//template<typename S>
//struct ShapeDistanceBisectionImpl<S, Box<S>, Sphere<S>>
//{
//  static bool run(
//      const BisectionSolver<S>& /*gjkSolver*/,
//      const Box<S>& s1,
//      const Transform3<S>& tf1,
//      const Sphere<S>& s2,
//      const Transform3<S>& tf2,
//      S* dist,
//      Vector3<S>* p1,
//      Vector3<S>* p2)
//  {
//    return detail::sphereBoxDistance(s2, tf2, s1, tf1, dist, p2, p1);
//  }
//};
//
////==============================================================================
//template<typename S>
//struct ShapeDistanceBisectionImpl<S, Sphere<S>, Capsule<S>>
//{
//  static bool run(
//      const BisectionSolver<S>& /*gjkSolver*/,
//      const Sphere<S>& s1,
//      const Transform3<S>& tf1,
//      const Capsule<S>& s2,
//      const Transform3<S>& tf2,
//      S* dist,
//      Vector3<S>* p1,
//      Vector3<S>* p2)
//  {
//    return detail::sphereCapsuleDistance(s1, tf1, s2, tf2, dist, p1, p2);
//  }
//};
//
////==============================================================================
//template<typename S>
//struct ShapeDistanceBisectionImpl<S, Capsule<S>, Sphere<S>>
//{
//  static bool run(
//      const BisectionSolver<S>& /*gjkSolver*/,
//      const Capsule<S>& s1,
//      const Transform3<S>& tf1,
//      const Sphere<S>& s2,
//      const Transform3<S>& tf2,
//      S* dist,
//      Vector3<S>* p1,
//      Vector3<S>* p2)
//  {
//    return detail::sphereCapsuleDistance(s2, tf2, s1, tf1, dist, p2, p1);
//  }
//};
//
////==============================================================================
//template<typename S>
//struct ShapeDistanceBisectionImpl<S, Sphere<S>, Cylinder<S>>
//{
//  static bool run(
//      const BisectionSolver<S>& /*gjkSolver*/,
//      const Sphere<S>& s1,
//      const Transform3<S>& tf1,
//      const Cylinder<S>& s2,
//      const Transform3<S>& tf2,
//      S* dist,
//      Vector3<S>* p1,
//      Vector3<S>* p2)
//  {
//    return detail::sphereCylinderDistance(s1, tf1, s2, tf2, dist, p1, p2);
//  }
//};
//
////==============================================================================
//template<typename S>
//struct ShapeDistanceBisectionImpl<S, Cylinder<S>, Sphere<S>>
//{
//  static bool run(
//      const BisectionSolver<S>& /*gjkSolver*/,
//      const Cylinder<S>& s1,
//      const Transform3<S>& tf1,
//      const Sphere<S>& s2,
//      const Transform3<S>& tf2,
//      S* dist,
//      Vector3<S>* p1,
//      Vector3<S>* p2)
//  {
//    return detail::sphereCylinderDistance(s2, tf2, s1, tf1, dist, p2, p1);
//  }
//};
//
////==============================================================================
//template<typename S>
//struct ShapeDistanceBisectionImpl<S, Sphere<S>, Sphere<S>>
//{
//  static bool run(
//      const BisectionSolver<S>& /*gjkSolver*/,
//      const Sphere<S>& s1,
//      const Transform3<S>& tf1,
//      const Sphere<S>& s2,
//      const Transform3<S>& tf2,
//      S* dist,
//      Vector3<S>* p1,
//      Vector3<S>* p2)
//  {
//    return detail::sphereSphereDistance(s1, tf1, s2, tf2, dist, p1, p2);
//  }
//};
//
////==============================================================================
//template<typename S>
//struct ShapeDistanceBisectionImpl<S, Capsule<S>, Capsule<S>>
//{
//  static bool run(
//      const BisectionSolver<S>& /*gjkSolver*/,
//      const Capsule<S>& s1,
//      const Transform3<S>& tf1,
//      const Capsule<S>& s2,
//      const Transform3<S>& tf2,
//      S* dist,
//      Vector3<S>* p1,
//      Vector3<S>* p2)
//  {
//    return detail::capsuleCapsuleDistance(s1, tf1, s2, tf2, dist, p1, p2);
//  }
//};

//==============================================================================
template<typename S, typename Shape>
struct ShapeTriangleDistanceBisectionImpl
{
  static bool run(
      const BisectionSolver<S>& gjkSolver,
      const Shape& s,
      const Transform3<S>& tf,
      const Vector3<S>& P1,
      const Vector3<S>& P2,
      const Vector3<S>& P3,
      S* dist,
      Vector3<S>* p1,
      Vector3<S>* p2)
  {
    TriangleP<S> tri(P1, P2, P3);
    detail::MinkowskiDiff<S> shape;
    shape.shapes[0] = &s;
    shape.shapes[1] = &tri;
    shape.toshape1 = tf.linear();
    shape.toshape0 = tf.inverse(Eigen::Isometry);

    bool res = detail::SeparationDistance(
          shape,
          gjkSolver.max_distance_iterations,
          gjkSolver.distance_tolerance,
          dist,
          p1,
          p2);

    return res;
  }
};

//==============================================================================
template<typename S>
template<typename Shape>
bool BisectionSolver<S>::shapeTriangleDistance(
    const Shape& s,
    const Transform3<S>& tf,
    const Vector3<S>& P1,
    const Vector3<S>& P2,
    const Vector3<S>& P3,
    S* dist,
    Vector3<S>* p1,
    Vector3<S>* p2) const
{
  return ShapeTriangleDistanceBisectionImpl<S, Shape>::run(
        *this, s, tf, P1, P2, P3, dist, p1, p2);
}

//==============================================================================
template<typename S>
struct ShapeTriangleDistanceBisectionImpl<S, Sphere<S>>
{
  static bool run(
      const BisectionSolver<S>& /*gjkSolver*/,
      const Sphere<S>& s,
      const Transform3<S>& tf,
      const Vector3<S>& P1,
      const Vector3<S>& P2,
      const Vector3<S>& P3,
      S* dist,
      Vector3<S>* p1,
      Vector3<S>* p2)
  {
    return detail::sphereTriangleDistance(s, tf, P1, P2, P3, dist, p1, p2);
  }
};

//==============================================================================
template<typename S, typename Shape>
struct ShapeTransformedTriangleDistanceBisectionImpl
{
  static bool run(
      const BisectionSolver<S>& gjkSolver,
      const Shape& s,
      const Transform3<S>& tf1,
      const Vector3<S>& P1,
      const Vector3<S>& P2,
      const Vector3<S>& P3,
      const Transform3<S>& tf2,
      S* dist,
      Vector3<S>* p1,
      Vector3<S>* p2)
  {
    TriangleP<S> tri(P1, P2, P3);

    detail::MinkowskiDiff<S> shape;
    shape.shapes[0] = &s;
    shape.shapes[1] = &tri;
    shape.toshape1.noalias() = tf2.linear().transpose() * tf1.linear();
    shape.toshape0 = tf1.inverse(Eigen::Isometry) * tf2;

    bool res = detail::SeparationDistance(
          shape,
          gjkSolver.max_distance_iterations,
          gjkSolver.distance_tolerance,
          dist,
          p1,
          p2);

    return res;
  }
};

//==============================================================================
template<typename S>
template<typename Shape>
bool BisectionSolver<S>::shapeTriangleDistance(
    const Shape& s,
    const Transform3<S>& tf1,
    const Vector3<S>& P1,
    const Vector3<S>& P2,
    const Vector3<S>& P3,
    const Transform3<S>& tf2,
    S* dist,
    Vector3<S>* p1,
    Vector3<S>* p2) const
{
  return ShapeTransformedTriangleDistanceBisectionImpl<S, Shape>::run(
        *this, s, tf1, P1, P2, P3, tf2, dist, p1, p2);
}

//==============================================================================
template<typename S>
struct ShapeTransformedTriangleDistanceBisectionImpl<S, Sphere<S>>
{
  static bool run(
      const BisectionSolver<S>& /*gjkSolver*/,
      const Sphere<S>& s,
      const Transform3<S>& tf1,
      const Vector3<S>& P1,
      const Vector3<S>& P2,
      const Vector3<S>& P3,
      const Transform3<S>& tf2,
      S* dist,
      Vector3<S>* p1,
      Vector3<S>* p2)
  {
    return detail::sphereTriangleDistance(
          s, tf1, P1, P2, P3, tf2, dist, p1, p2);
  }
};

//==============================================================================
//==============================================================================
template<typename S, typename Shape1, typename Shape2>
struct ShapeSignedDistanceBisectionImpl
{
  static bool run(
      const BisectionSolver<S>& gjkSolver,
      const Shape1& s1,
      const Transform3<S>& tf1,
      const Shape2& s2,
      const Transform3<S>& tf2,
      S* dist,
      Vector3<S>* p1,
      Vector3<S>* p2)
  {
    detail::MinkowskiDiff<S> shape;
    shape.shapes[0] = &s1;
    shape.shapes[1] = &s2;
    shape.toshape1.noalias() = tf2.linear().transpose() * tf1.linear();
    shape.toshape0 = tf1.inverse(Eigen::Isometry) * tf2;

    bool res = detail::SignedDistance(
          shape,
          gjkSolver.max_distance_iterations,
          gjkSolver.distance_tolerance,
          dist,
          p1,
          p2);

    return res;
  }
};

template<typename S>
template<typename Shape1, typename Shape2>
bool BisectionSolver<S>::shapeSignedDistance(
    const Shape1& s1,
    const Transform3<S>& tf1,
    const Shape2& s2,
    const Transform3<S>& tf2,
    S* dist,
    Vector3<S>* p1,
    Vector3<S>* p2) const
{
  bool result = false;
  try {
    result = ShapeSignedDistanceBisectionImpl<S, Shape1, Shape2>::run(
        *this, s1, tf1, s2, tf2, dist, p1, p2);
  } catch (const FailedAtThisConfiguration& e) {
    ThrowDetailedConfiguration(s1, tf1, s2, tf2, *this, e);
  }
  return result;
}

//==============================================================================
template<typename S>
BisectionSolver<S>::BisectionSolver()
{
  max_collision_iterations = 500;
  max_distance_iterations = 1000;
  collision_tolerance = constants<S>::gjk_default_tolerance();
  distance_tolerance = 1e-6;
}

//==============================================================================
template<typename S>
void BisectionSolver<S>::enableCachedGuess(bool if_enable) const
{
  FCL_UNUSED(if_enable);

  // TODO: need change libccd to exploit spatial coherence
}

//==============================================================================
template<typename S>
void BisectionSolver<S>::setCachedGuess(
    const Vector3<S>& guess) const
{
  FCL_UNUSED(guess);

  // TODO: need change libccd to exploit spatial coherence
}

//==============================================================================
template<typename S>
Vector3<S> BisectionSolver<S>::getCachedGuess() const
{
  return Vector3<S>(-1, 0, 0);
}

} // namespace detail
} // namespace fcl

#endif
