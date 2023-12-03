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

#ifndef FCL_NARROWPHASE_DETAIL_BISECTIONDISTANCE_INL_H
#define FCL_NARROWPHASE_DETAIL_BISECTIONDISTANCE_INL_H

#include "fcl/narrowphase/detail/convexity_based_algorithm/bisection_distance.h"
#include "fcl/narrowphase/detail/failed_at_this_configuration.h"

#include <array>
#include <unordered_map>
#include <unordered_set>

#include "fcl/common/unused.h"
#include "fcl/common/warning.h"

namespace fcl
{

namespace detail
{

//==============================================================================
extern template
bool GJKCollide(const MinkowskiDiffd& shape, unsigned int max_iterations, double tolerance,
    Vector3d* contact_points, double* penetration_depth, Vector3d* normal);

//==============================================================================
extern template
bool SeparationDistance(const MinkowskiDiffd& shape, unsigned int max_iterations, double tolerance,
    double* dist, Vector3d* p1, Vector3d* p2);

extern template
bool SignedDistance(const MinkowskiDiffd& shape, unsigned int max_iterations, double tolerance,
    double* dist, Vector3d* p1, Vector3d* p2);

namespace libccd_extension
{

// Compares the given `value` against a _squared epsilon_. This is particularly
// important when testing some quantity (e.g., distance) to see if it
// is _functionally_ zero but using its _squared_ value in the test. Comparing
// _squared distance_ directly against epsilon is equivalent to comparing
// distance to sqrt(epsilon) -- we classify the distance as zero or not using
// only half the available precision.
template <typename S>
static bool isAbsValueLessThanEpsSquared(S val)
{
    return std::abs(val) < std::numeric_limits<S>::epsilon() * std::numeric_limits<S>::epsilon();
}

/** Determines if the the triangle defined by the three vertices has zero area.
 Area can be zero for one of two reasons:
   - the triangle is so small that the vertices are functionally coincident, or
   - the vertices are co-linear.
 Both conditions are computed with respect to machine precision.
 @returns true if the area is zero.  */
template <typename S>
static bool triangle_area_is_zero(const Vector3<S>& a, const Vector3<S>& b, const Vector3<S>& c)
{
  // First coincidence condition. This doesn't *explicitly* test for b and c
  // being coincident. That will be captured in the subsequent co-linearity
  // test. If b and c *were* coincident, it would be cheaper to perform the
  // coincidence test than the co-linearity test.
  // However, the expectation is that typically the triangle will not have zero
  // area. In that case, we want to minimize the number of tests performed on
  // the average, so we prefer to eliminate one coincidence test.
  if (a.isApprox(b) || a.isApprox(c))
      return true;

  // We're going to compute the *sine* of the angle θ between edges (given that
  // the vertices are *not* coincident). If the sin(θ) < ε, the edges are
  // co-linear.
  Vector3<S> AB = (b - a).normalized(), AC = (c - a).normalized(), n = AB.cross(AC);
  constexpr S eps = constants<S>::eps();
  // Second co-linearity condition.
  if (n.squaredNorm() < eps * eps)
      return true;
  return false;
}

template <typename S>
static S simplexReduceToTriangle(Simplex<S> &simplex, S &dist, Vector3<S> &best_witness)
{
  S newdist;
  Vector3<S> witness;
  int best = -1;

  // try the fourth point in all three positions
  for (int i = 0; i < 3; i++)
  {
    newdist = Project<S>::originTriDist2(simplex.at(i == 0 ? 3 : 0).v, simplex.at(i == 1 ? 3 : 1).v, simplex.at(i == 2 ? 3 : 2).v, witness);
    newdist = std::sqrt(newdist);

    // record the best triangle
    if (newdist < dist)
    {
      dist = newdist;
      best = i;
      best_witness = witness;
    }
  }

  if (best >= 0)
    simplex.set(best, simplex.at(3));
  simplex.setSize(3);

  return dist;
}

template <typename S>
inline void tripleCross(const Vector3<S> &a, const Vector3<S> &b, const Vector3<S> &c, Vector3<S> &d)
{
  d = a.cross(b).cross(c);
}

/* This is *not* an implementation of the general function: what's the nearest
 point on the line segment AB to the origin O? It is not intended to be.
 This is a limited, special case which exploits the known (or at least
 expected) construction history of AB. The history is as follows:

   1. We originally started with the Minkowski support point B (p_OB), which
      was *not* the origin.
   2. We define a support direction as p_BO = -p_OB and use that to get the
      Minkowski support point A.
   3. We confirm that O is not strictly beyond A in the direction p_BO
      (confirming separation).
   4. Then, and only then, do we invoke this method.

 The method will do one of two things:

   - determine if the origin lies within the simplex (i.e. lies on the line
     segment, confirming non-separation) and reports if this is the case,
   - otherwise it computes a new support direction: a vector pointing to the
     origin from the nearest point on the segment AB. The direction is
     guaranteed; the only guarantee about the magnitude is that it is
     numerically viable (i.e. greater than epsilon).

 The algorithm exploits the construction history as outlined below. Without
 loss of generality, we place B some non-zero distance away from the origin
 along the î direction (all other orientations can be rigidly transformed to
 this canonical example). The diagram below shows the origin O and the point
 B. It also shows three regions: 1, 2, and 3.

                     ĵ
              1      ⯅    2        3
                     │░░░░░░░░░░┆▒▒▒▒▒
                     │░░░░░░░░░░┆▒▒▒▒▒
                     │░░░░░░░░░░┆▒▒▒▒▒
                     │░░░░░░░░░░┆▒▒▒▒▒
                     │░░░░░░░░░░┆▒▒▒▒▒
          ───────────O──────────B────⯈  î
                     │░░░░░░░░░░┆▒▒▒▒▒
                     │░░░░░░░░░░┆▒▒▒▒▒
                     │░░░░░░░░░░┆▒▒▒▒▒
                     │░░░░░░░░░░┆▒▒▒▒▒
                     │░░░░░░░░░░┆▒▒▒▒▒

 The point A cannot lie in region 3.

   - B is a support point of the Minkowski difference. p_BO defines the
     support vector that produces the support point A. Vertex A must lie at
     least as far in the p_BO as B otherwise A is not actually a valid support
     point for that direction. It could lie on the boundary of regions 2 & 3
     and still be a valid support point.

 The point A cannot lie in region 2 (or on the boundary between 2 & 3).
   - We confirm that O is not strictly beyond A in the direction p_BO. For all
     A in region 2, O lies beyond A (when projected onto the p_BO vector).

 The point A _must_ lie in region 1 (including the boundary between regions 1 &
 2) by process of elimination.

 The implication of this is that the O must project onto the _interior_ of the
 line segment AB (with one notable exception). If A = O, then the projection of
 O is at the end point and, is in fact, itself.

 Therefore, this function can only have two possible outcomes:

   1. In the case where p_BA = k⋅p_BO (i.e., they are co-linear), we know the
      origin lies "in" the simplex. If A = O, it lies on the simplex's surface
      and the objects are touching, otherwise, the objects are penetrating.
      Either way, we can report that they are definitely *not* separated.
   2. p_BA ≠ k⋅p_BO, we define the new support direction as perpendicular to the
      line segment AB, pointing to O from the nearest point on the segment to O.

 Return value indicates concrete knowledge that the origin lies "in" the
 2-simplex (encoded as a 1), or indication that computation should continue (0).

 @note: doSimplex2 should _only_ be called with the construction history
 outlined above: promotion of a 1-simplex. doSimplex2() is only invoked by
 doSimplex(). This follows the computation of A and the promotion of the
 simplex. Therefore, the history is always valid. Even though doSimplex3() can
 demote itself to a 2-simplex, that 2-simplex immediately gets promoted back to
 a 3-simplex via the same construction process. Therefore, as long as
 doSimplex2() is only called from doSimplex() its region 1 assumption _should_
 remain valid.
*/
template <typename S>
static int doSimplex2(const Support<S> &a, const Support<S> &b, Vector3<S> &dir)
{
  // Used to define numerical thresholds near zero; typically scaled to the size
  // of the quantities being tested.
  constexpr S eps = constants<S>::eps();

  const Vector3<S> p_OA(a.v);
  const Vector3<S> p_OB(b.v);

  // Confirm that A is in region 1. Given that A may be very near to the origin,
  // we must avoid normalizing p_OA. So, we use this instead.
  //  let A' be the projection of A onto the line defined by O and B.
  //  |A'B| >= |OB| iff A is in region 1.
  // Numerically, we can express it as follows (allowing for |A'B| to be ever
  // so slightly *less* than |OB|):
  //             p_AB ⋅ phat_OB >= |p_OB| - |p_OB| * ε = |p_OB| * (1 - ε)
  //    p_AB ⋅ phat_OB ⋅ |p_OB| >= |p_OB|⋅|p_OB| * (1 - ε)
  //                p_AB ⋅ p_OB >= |p_OB|² * (1 - ε)
  //       (p_OB - p_OA) ⋅ p_OB >= |p_OB|² * (1 - ε)
  //  p_OB ⋅ p_OB - p_OA ⋅ p_OB >= |p_OB|² * (1 - ε)
  //      |p_OB|² - p_OA ⋅ p_OB >= |p_OB|² * (1 - ε)
  //               -p_OA ⋅ p_OB >= -|p_OB|²ε
  //                p_OA ⋅ p_OB <= |p_OB|²ε
  assert(p_OA.dot(p_OB) <= p_OB.squaredNorm() * eps && "A is not in region 1");

  // Test for co-linearity. Given A is in region 1, co-linearity --> O is "in"
  // the simplex.
  // We'll use the angle between two vectors to determine co-linearity: p_AB
  // and p_OB. If they are co-linear, then the angle between them (θ) is zero.
  // Similarly, sin(θ) is zero. Ideally, it can be expressed as:
  //   |p_AB × p_OB| = |p_AB||p_OB||sin(θ)| = 0
  // Numerically, we allow θ (and sin(θ)) to be slightly larger than zero
  // leading to a numerical formulation as:
  //   |p_AB × p_OB| = |p_AB||p_OB||sin(θ)| < |p_AB||p_OB|ε
  // Finally, to reduce the computational cost, we eliminate the square roots by
  // evaluating the equivalently discriminating test:
  //   |p_AB × p_OB|² < |p_AB|²|p_OB|²ε²
  //
  // In addition to providing a measure of co-linearity, the cross product gives
  // us the normal to the plane on which points A, B, and O lie (which we will
  // use later to compute a new support direction, as necessary).
  const Vector3<S> p_AB = p_OB - p_OA;
  const Vector3<S> plane_normal = p_OB.cross(p_AB);
  if (plane_normal.squaredNorm() < p_AB.squaredNorm() * p_OB.squaredNorm() * eps * eps)
    return 1;

  // O is not co-linear with AB, so dist(O, AB) > ε. Define `dir` as the
  // direction to O from the nearest point on AB.
  // Note: We use the normalized `plane_normal` (n̂) because we've already
  // concluded that the origin is farther from AB than ε. We want to make sure
  // `dir` likewise has a magnitude larger than ε. With normalization, we know
  // |dir| = |n̂ × AB| = |AB| > dist(O, AB) > ε.
  // Without normalizing, if |OA| and |OB| were smaller than ³√ε but
  // sufficiently larger than ε, dist(O, AB) > ε, but |dir| < ε.
  dir = plane_normal.normalized().cross(p_AB);
  return 0;
}

template <typename S>
static int doSimplex2(Simplex<S> &simplex, Vector3<S> &dir)
{
  return doSimplex2(simplex.at(1), simplex.at(0), dir);
}

// TODO(SeanCurtis-TRI): Define the return value:
//   1: (like doSimplex2) --> origin is "in" the simplex.
//   0:
//  -1: If the 3-simplex is degenerate. How is this intepreted?
template <typename S>
static int doSimplex3(Simplex<S> &simplex, Vector3<S> &dir)
{
  const Support<S> &A = simplex.at(2), &B = simplex.at(1), &C = simplex.at(0);

  // check touching contact
  // Compute origin_projection as well. Without computing the origin projection,
  // libccd could give inaccurate result. See
  // https://github.com/danfis/libccd/issues/55.
  Vector3<S> origin_projection_unused;
  const S dist_squared = Project<S>::originTriDist2(A.v, B.v, C.v, origin_projection_unused);
  if (isAbsValueLessThanEpsSquared(dist_squared))
    return 1;

  // check if triangle is really triangle (has area > 0)
  // if not simplex can't be expanded and thus no intersection is found
  // TODO(SeanCurtis-TRI): Coincident points is sufficient but not necessary
  // for a zero-area triangle. What about co-linearity? Can we guarantee that
  // co-linearity can't happen?  See the `triangle_area_is_zero()` method in
  // this same file.
  if (A.v.isApprox(B.v) || A.v.isApprox(C.v))
    // TODO(SeanCurtis-TRI): Why do we simply quit if the simplex is degenerate?
    return -1;

  Vector3<S> AO = -A.v, AB = B.v - A.v, AC = C.v - A.v;
  Vector3<S> ABC = AB.cross(AC), tmp = ABC.cross(AC);
  S dot = tmp.dot(AO), eps = std::numeric_limits<S>::epsilon();

  if (abs(dot) < eps || dot > 0.0)
  {
    dot = AC.dot(AO);
    if (abs(dot) < eps || dot > 0.0)
    {
      // C is already in place
      simplex.set(1, A);
      simplex.setSize(2);
      tripleCross(AC, AO, AC, dir);
    }
    else
    {
    ccd_do_simplex3_45:
      dot = AB.dot(AO);
      if (abs(dot) < eps || dot > 0.0)
      {
        simplex.set(0, B);
        simplex.set(1, A);
        simplex.setSize(2);
        tripleCross(AB, AO, AB, dir);
      }
      else
      {
        simplex.set(0, A);
        simplex.setSize(1);
        dir = AO;
      }
    }
  }
  else
  {
    tmp = AB.cross(ABC);
    dot = tmp.dot(AO);
    if (abs(dot) < eps || dot > 0.0)
    {
      goto ccd_do_simplex3_45;
    }
    else
    {
      dot = ABC.dot(AO);
      if (abs(dot) < eps || dot > 0.0)
        dir = ABC;
      else
      {
        Support<S> Ctmp = C;
        simplex.set(0, B);
        simplex.set(1, Ctmp);
        dir = -ABC;
      }
    }
  }

  return 0;
}

template <typename S>
static int doSimplex4(Simplex<S> &simplex, Vector3<S> &dir)
{
  const Support<S> &A = simplex.at(3), &B = simplex.at(2), &C = simplex.at(1), &D = simplex.at(0);

  // check if tetrahedron is really tetrahedron (has volume > 0)
  // if it is not simplex can't be expanded and thus no intersection is
  // found.
  // point_projection_on_triangle_unused is not used. We ask
  // ccdVec3PointTriDist2 to compute this witness point, so as to get a
  // numerical robust dist_squared. See
  // https://github.com/danfis/libccd/issues/55 for an explanation.
  Vector3<S> point_projection_on_triangle_unused;
  S dist_squared = Project<S>::pointTriDist2(B.v, C.v, D.v, A.v, point_projection_on_triangle_unused);
  if (isAbsValueLessThanEpsSquared(dist_squared))
    return -1;

  // check if origin lies on some of tetrahedron's face - if so objects
  // intersect
  dist_squared = Project<S>::originTriDist2(A.v, B.v, C.v, point_projection_on_triangle_unused);
  if (isAbsValueLessThanEpsSquared((dist_squared))) return 1;
  dist_squared = Project<S>::originTriDist2(A.v, C.v, D.v, point_projection_on_triangle_unused);
  if (isAbsValueLessThanEpsSquared((dist_squared))) return 1;
  dist_squared = Project<S>::originTriDist2(A.v, B.v, D.v, point_projection_on_triangle_unused);
  if (isAbsValueLessThanEpsSquared((dist_squared))) return 1;
  dist_squared = Project<S>::originTriDist2(B.v, C.v, D.v, point_projection_on_triangle_unused);
  if (isAbsValueLessThanEpsSquared((dist_squared))) return 1;

  Vector3<S> AO = -A.v, AB = B.v - A.v, AC = C.v - A.v, AD = D.v - A.v;
  Vector3<S> ABC= AB.cross(AC), ACD = AC.cross(AD), ADB = AD.cross(AB);

  auto ccdSign = [](S x) -> int {
    return abs(x) < std::numeric_limits<S>::epsilon() ? 0 : (x > 0 ? 1 : -1);
  };

  // side (positive or negative) of B, C, D relative to planes ACD, ADB
  // and ABC respectively
  int B_on_ACD = ccdSign(ACD.dot(AB));
  int C_on_ADB = ccdSign(ADB.dot(AC));
  int D_on_ABC = ccdSign(ABC.dot(AD));

  // whether origin is on same side of ACD, ADB, ABC as B, C, D
  // respectively
  int AB_O = ccdSign(ACD.dot(AO)) == B_on_ACD;
  int AC_O = ccdSign(ADB.dot(AO)) == C_on_ADB;
  int AD_O = ccdSign(ABC.dot(AO)) == D_on_ABC;

  if (AB_O && AC_O && AD_O)
    // origin is in tetrahedron
    return 1;
    // rearrange simplex to triangle and call doSimplex3()
  else if (!AB_O)
  {
    // B is farthest from the origin among all of the tetrahedron's
    // points, so remove it from the list and go on with the triangle
    // case

    // D and C are in place
    simplex.set(2, A);
    simplex.setSize(3);
  }
  else if (!AC_O)
  {
    // C is farthest
    simplex.set(1, D);
    simplex.set(0, B);
    simplex.set(2, A);
    simplex.setSize(3);
  }
  else
  { // (!AD_O)
    simplex.set(0, C);
    simplex.set(1, B);
    simplex.set(2, A);
    simplex.setSize(3);
  }

  return doSimplex3(simplex, dir);
}

template <typename S>
static int doSimplex(Simplex<S> &simplex, Vector3<S> &dir)
{
  if (simplex.size() == 2)
    // simplex contains segment only one segment
    return doSimplex2(simplex, dir);
  else if (simplex.size() == 3)
    // simplex contains triangle
    return doSimplex3(simplex, dir);
  else // ccdSimplexSize(simplex) == 4
    // tetrahedron - this is the only shape which can encapsule origin
    // so doSimplex4() also contains test on it
    return doSimplex4(simplex, dir);
}

template <typename S>
static int __ccdGJK(const MinkowskiDiff<S>& shape, Simplex<S> &simplex, unsigned int max_iterations)
{
  Vector3<S> dir(1.0, 0.0, 0.0); // direction vector
  Support<S> last = shape.supportS(dir); // last support point
  simplex.add(last);
  if (last.v.dot(dir) < 0.0)
    return -1; // intersection not found
  // set up direction vector to as (O - last) which is exactly -last
  dir = (-last.v).normalized();
  constexpr S eps = constants<S>::eps();
  // start iterations
  for (unsigned int iterations = 0; iterations < max_iterations; ++iterations)
  {
    // obtain support point
    last = shape.supportS(dir);

    // add last support vector to simplex
    simplex.add(last);

    // check if farthest point in Minkowski difference in direction dir
    // isn't somewhere before origin (the test on negative dot product)
    // - because if it is, objects are not intersecting at all.
    if (last.v.dot(dir) < 0.0)
      return -1; // intersection not found

    // if doSimplex returns 1 if objects intersect, -1 if objects don't
    // intersect and 0 if algorithm should continue
    int do_simplex_res = doSimplex(simplex, dir);
    if (do_simplex_res == 1)
      return 0; // intersection found
    if (do_simplex_res == -1)
      return -1; // intersection not found
    if (dir.squaredNorm() < eps)
      return -1; // intersection not found
    dir = dir.normalized();
  }
  // intersection wasn't found
  return -1;
}

/** Given a single support point, `q`, extract the point `p1` and `p2`, the
 points on object 1 and 2, respectively, in the support data of `q`.  */
template <typename S>
static void extractObjectPointsFromPoint(const Support<S> &q, Vector3<S> &p1, Vector3<S> &p2)
{
  // TODO(SeanCurtis-TRI): Determine if I should be demanding that p1 and p2
  // are defined.
  // Closest points are the ones stored in the simplex
  p1 = q.v1;
  p2 = q.v2;
}

/** Given two support points which define a line segment (`a` and `b`), and a
 point on that line segment `p`, computes the points `p1` and `p2`, the points
 on object 1 and 2, respectively, in the support data which correspond to `p`.
 @pre `p = a + s(b - a), 0 <= s <= 1`  */
template <typename S>
static S alphaRatio(const Support<S> &a, const Support<S> &b, const Vector3<S> &p)
{
  // Closest points lie on the segment defined by the points in the simplex
  // Let the segment be defined by points A and B. We can write p as
  //
  // p = A + s*AB, 0 <= s <= 1
  // p - A = s*AB
  Vector3<S> AB = b.v - a.v;

  // This defines three equations, but we only need one. Taking the i-th
  // component gives
  //
  // p_i - A_i = s*AB_i.
  //
  // Thus, s is given by
  //
  // s = (p_i - A_i)/AB_i.
  //
  // To avoid dividing by an AB_i ≪ 1, we choose i such that |AB_i| is
  // maximized
  S abs_AB_x{std::abs(AB(0))};
  S abs_AB_y{std::abs(AB(1))};
  S abs_AB_z{std::abs(AB(2))};

  S A_i, AB_i, p_i;
  if (abs_AB_x >= abs_AB_y && abs_AB_x >= abs_AB_z)
  {
    A_i = a.v[0];
    AB_i = AB[0];
    p_i = p[0];
  }
  else if (abs_AB_y >= abs_AB_z)
  {
    A_i = a.v[1];
    AB_i = AB[1];
    p_i = p[1];
  }
  else
  {
    A_i = a.v[2];
    AB_i = AB[2];
    p_i = p[2];
  }

  if (std::abs(AB_i) < constants<S>::eps())
      return -1.0;

  // TODO(SeanCurtis-TRI): If p1 or p2 is null, there seems little point in
  // calling this method. It seems that both of these being non-null should be
  // a *requirement*. Determine that this is the case and do so.
  S s = (p_i - A_i) / AB_i;
  return s;
}

template <typename S>
static void extractObjectPointsFromSegment(const Support<S> &a, const Support<S> &b,
                                           Vector3<S> &p1, Vector3<S> &p2,
                                           const Vector3<S> &p)
{
  S s = alphaRatio(a, b, p);
  if (s < 0.0) 
  {
    // Points are coincident; treat as a single point.
    extractObjectPointsFromPoint(a, p1, p2);
    return;
  }

  auto calc_p = [](const Vector3<S> &p_a, const Vector3<S> &p_b, Vector3<S> &p, S s) {
    p = p_a + s * (p_b - p_a);
  };

  calc_p(a.v1, b.v1, p1, s);
  calc_p(a.v2, b.v2, p2, s);
}

/** Returns the points `p1` and `p2` on the original shapes that correspond to
 point `p` in the given simplex.
 @pre simplex size <= 3.
 @pre p lies _on_ the simplex (i.e., within the triangle, line segment, or is
      coincident with the point).  */
template <typename S>
static void extractClosestPoints(const Simplex<S> &simplex, Vector3<S> &p1, Vector3<S> &p2, const Vector3<S> &p)
{
  const int simplex_size = simplex.size();
  assert(simplex_size <= 3);
  if (simplex_size == 1)
    extractObjectPointsFromPoint(simplex.at(0), p1, p2);
  else if (simplex_size == 2)
    extractObjectPointsFromSegment(simplex.at(0), simplex.at(1), p1, p2, p);
  else // simplex_size == 3
  {
    if (triangle_area_is_zero(simplex.at(0).v, simplex.at(1).v, simplex.at(2).v))
    {
      // The triangle is degenerate; compute the nearest point to a line
      // segment. The segment is defined by the most-distant vertex pair.
      int a_index, b_index;
      S AB_len2 = (simplex.at(1).v - simplex.at(0).v).squaredNorm();
      S AC_len2 = (simplex.at(2).v - simplex.at(0).v).squaredNorm();
      S BC_len2 = (simplex.at(2).v - simplex.at(1).v).squaredNorm();

      if (AB_len2 >= AC_len2 && AB_len2 >= BC_len2)
      {
        a_index = 0;
        b_index = 1;
      }
      else if (AC_len2 >= AB_len2 && AC_len2 >= BC_len2)
      {
        a_index = 0;
        b_index = 2;
      }
      else
      {
        a_index = 1;
        b_index = 2;
      }
      extractObjectPointsFromSegment(simplex.at(a_index), simplex.at(b_index), p1, p2, p);
      return;
    }

    // Compute the barycentric coordinates of point p in triangle ABC.
    //
    //             A
    //             ╱╲                p = αA + βB + γC
    //            ╱ |╲
    //           ╱  | ╲              α = 1 - β - γ
    //          ╱ p |  ╲             β = AREA(pAC) / AREA(ABC)
    //         ╱   / \  ╲            γ = AREA(pAB) / AREA(ABC)
    //        ╱__/_____\_╲
    //      B             C          AREA(XYZ) = |r_XY × r_XZ| / 2
    //
    //  Rewrite coordinates in terms of cross products.
    //
    //    β = AREA(pAC) / AREA(ABC) = |r_Ap × r_AC| / |r_AB × r_AC|
    //    γ = AREA(pAB) / AREA(ABC) = |r_AB × r_Ap| / |r_AB × r_AC|
    //
    // NOTE: There are multiple options for the cross products, these have been
    // selected to re-use as many symbols as possible.
    //
    // Solving for β and γ:
    //
    //  β = |r_Ap × r_AC| / |r_AB × r_AC|
    //  β = |r_Ap × r_AC| / |n|                  n ≙ r_AB × r_AC, n̂ = n / |n|
    //  β = n̂·(r_Ap × r_AC) / n̂·n                This step arises from the fact
    //                                           that (r_Ap × r_AC) and n point
    //                                           in the same direction. It
    //                                           allows us to take a single sqrt
    //                                           instead of three.
    //  β = (n/|n|)·(r_Ap × r_AC) / (n/|n|)·n
    //  β = n·(r_Ap × r_AC) / n·n
    //  β = n·(r_Ap × r_AC) / |n|²
    //
    // A similar process to solve for gamma
    //  γ = n·(r_AB × r_Ap) / |n|²

    // Compute n and |n|².
    Vector3<S> r_AB = simplex.at(1).v - simplex.at(0).v;
    Vector3<S> r_AC = simplex.at(2).v - simplex.at(0).v, n = r_AB.cross(r_AC);
    S norm_squared_n = n.squaredNorm();

    // Compute r_Ap.
    Vector3<S> r_Ap = p - simplex.at(0).v;

    // Compute the cross products in the numerators.
    Vector3<S> r_Ap_cross_r_AC = r_Ap.cross(r_AC);
    Vector3<S> r_AB_cross_r_Ap = r_AB.cross(r_Ap);

    // Compute beta and gamma.
    S beta = n.dot(r_Ap_cross_r_AC) / norm_squared_n;
    S gamma= n.dot(r_AB_cross_r_Ap) / norm_squared_n;

    // Evaluate barycentric interpolation (with the locally defined barycentric
    // coordinates).
    auto interpolate = [&beta, &gamma](const Vector3<S>& r_WA,
                                       const Vector3<S>& r_WB,
                                       const Vector3<S>& r_WC,
                                       Vector3<S> &r_WP)
    {
      // r_WP = r_WA + β * r_AB + γ * r_AC
      r_WP = r_WA + beta * (r_WB - r_WA) + gamma * (r_WC - r_WA);
    };

    interpolate(simplex.at(0).v1, simplex.at(1).v1, simplex.at(2).v1, p1);
    interpolate(simplex.at(0).v2, simplex.at(1).v2, simplex.at(2).v2, p2);
  }
}

template <typename S>
static void _simplexClosestP(Simplex<S> &simplex, Vector3<S> &closest_p, S &distp)
{
    int sz = simplex.size();
    if (sz == 1)
    {
        closest_p = simplex.at(0).v;
        distp = closest_p.norm();
    }
    else if (sz == 2)
    {
        distp = Project<S>::originSegmentDist2(simplex.at(0).v, simplex.at(1).v, closest_p);
        distp = std::sqrt(distp);
    }
    else if (sz == 3)
    {
        distp = Project<S>::originTriDist2(simplex.at(0).v, simplex.at(1).v, simplex.at(2).v, closest_p);
        distp = std::sqrt(distp);
    }
    else
        distp = simplexReduceToTriangle(simplex, distp, closest_p);
}

template <typename S>
static void _ccdSupport(const MinkowskiDiff<S>& shape, const Vector3<S> &closest_p,
        Vector3<S> &dir, Support<S> &last)
{
    // point direction towards the origin
    dir = (-closest_p).normalized();
    last= shape.supportS(dir);
}

template <typename S>
static bool _ccdOptimal(const Vector3<S> &closest_p, const Vector3<S> &dir, const Vector3<S> &last, S tol)
{
    if (dir.dot(last - closest_p) < tol)
        return true;
    return false;
}

// Computes the distance between two non-penetrating convex objects, returning the distance and
// nearest points on each object.
// @param simplex A witness to the objects' separation generated by the GJK
// algorithm. NOTE: the simplex is not necessarily sufficiently refined to
// report the actual distance and may be further refined in this method.
// @param p1 If the objects are non-penetrating, the point on the surface of
// obj1 closest to obj2 (expressed in the world frame).
// @param p2 If the objects are non-penetrating, the point on the surface of
// obj2 closest to obj1 (expressed in the world frame).
// @returns The minimum distance between the two objects. If they are
// penetrating, -1 is returned.
template <typename S>
static inline S _ccdDist(const MinkowskiDiff<S>& shape, Simplex<S> &simplex, unsigned int max_iterations, S tol,
        Vector3<S> &p1, Vector3<S> &p2)
{
    S last_dist = std::numeric_limits<S>::max();
    for (unsigned int iterations = 0u; iterations < max_iterations; ++iterations)
    {
        Vector3<S> closest_p; // The point on the simplex that is closest to the
        // origin.
        // get a next direction vector
        // we are trying to find out a point on the minkowski difference
        // that is nearest to the origin, so we obtain a point on the
        // simplex that is nearest and try to exapand the simplex towards
        // the origin
        _simplexClosestP(simplex, closest_p, last_dist);

        Vector3<S> dir; // direction vector
        Support<S> last; // last support point
        _ccdSupport(shape, closest_p, dir, last);

        if (_ccdOptimal(closest_p, dir, last.v, tol))
        {
            extractClosestPoints(simplex, p1, p2, closest_p);
            return last_dist;
        }

        // add a point to simplex
        simplex.add(last);
    }

    return -S(1.0);
}

template <typename S>
static inline S _ccdDistS(const MinkowskiDiff<S>& shape, Simplex<S> &simplex, unsigned int max_iterations, S tol,
        Vector3<S> &p1, Vector3<S> &p2)
{
    S last_dist = -std::numeric_limits<S>::max(), dist = last_dist, distp = -last_dist;
    Vector3<S> dir, best_dir; // direction vector
    Support<S> last, best_support; // last support point

    bool simplified = false;
    for (unsigned int iterations = 0u; iterations < max_iterations; ++iterations)
    {
        if (simplified)
            _ccdSupport(shape, best_support.v, dir, last);
        else 
        {
            Vector3<S> closest_p; 
            _simplexClosestP(simplex, closest_p, distp);
            _ccdSupport(shape, closest_p, dir, last);
            if (_ccdOptimal(closest_p, dir, last.v, tol)) // termination one 
            {
                extractClosestPoints(simplex, p1, p2, closest_p);
                return distp;
            }
        }

        dist = -dir.dot(last.v);
        bool invalid = std::max(dist, last_dist) <= 0.0;
        if (invalid || iterations == 0u)
        {
            simplified = true;
            if (dist > last_dist)
            {
                last_dist = dist;
                best_dir = dir;
                best_support = last;
            }
            if (invalid)
            {
                simplified = false;
                simplex.add(last);
            }
            continue;
        }

        Vector3<S> p_B, p_BA;
        bool right = true, bisection = false;
        if (dist > last_dist)
        {
            p_B = dir; p_BA = best_dir - p_B;
            bisection = p_BA.dot(last.v) < -tol;
        }
        else
        {
            right = false;
            p_B = best_dir; p_BA = dir - p_B;
            bisection = p_BA.dot(best_support.v) < -tol;
        }

        if (bisection) // core part
        {
            S last_dist = right ? dist : last_dist;
            Vector3<S> dir1 = (p_B + tol * p_BA).normalized();
            Support<S> temp = shape.supportS(dir1);
            S d = -dir1.dot(temp.v);
            if (d > last_dist)
            {
                dist = last_dist = d;
                dir = dir1; last = temp;

                S s = 0.5, d1 = -std::numeric_limits<S>::max();
                while (tol < s)
                {
                    dir1 = (p_B + s * p_BA).normalized();
                    temp = shape.supportS(dir1);
                    d = -dir1.dot(temp.v);
                    if (d > last_dist)
                    {
                        dist = last_dist = d1 = d;
                        dir = dir1; last = temp; s *= 0.5;
                        continue;
                    }
                    if (d <= d1)
                        break;
                    else 
                        d1 = d; 
                    //if (temp.v.isApprox(best_support.v) || temp.v.isApprox(last.v)) // box
                    //    break;
                    s *= 0.5;
                }
            }
        }

        if (dist >= last_dist && dist - last_dist < tol) // termination two // todo
        {
            p1 = shape.support0(dir).dot(dir)*dir;
            p2 = shape.support1(-dir).dot(dir)*dir;
            return dist;
        }
        if (dist > last_dist)
        {
            simplified = true;
            best_dir = dir;
            best_support = last;
            last_dist = dist;
        }
        else
        {
            if (simplified)
            {
                simplex.set(0, best_support);
                simplex.setSize(1);
            }
            simplex.add(last);
            simplified = false;
        }
    }

    return last_dist;
}

template <typename S>
static inline bool _bisectionOptimal(const Vector3<S> &dir, const Vector3<S> &point, S tol, Vector3<S> &normal)
{
    normal = point.cross(dir);
    S len = normal.norm();
    if (len < tol)
    {
        //std::cout << "Optimal found!" << std::endl;
        return true;
    }
    normal /= len;
    return false;
}

// A negative return value means separation
template <typename S>
static inline Vector3<S> _bisectionDistanceCore(const MinkowskiDiff<S> & shape, unsigned int max_iterations, S tol, unsigned int & iterations,
    S & dist_best, bool last_best, Vector3<S> & last, Vector3<S> & current, const Vector3<S> & normal, Vector3<S> & p1, Vector3<S> & p2)
{
    Vector3<S> dir_best;
    for (; iterations < max_iterations; ++iterations) // bisection core
    {
        Vector3<S> dir = normal.cross(last - current).normalized();
        Vector3<S> temp = shape.support(dir);
        S dist = dir.dot(temp);
        if (dist > dist_best)
        {
            if (last_best)
                current = temp;
            else
                last = temp;
            continue;
        }
        dir_best = dir;
        dist_best = dist;
        if (dist - dir.dot(last) < tol) break;
        Vector3<S> tnormal = temp.cross(dir);
        if (tnormal.dot(normal) > S(0.0))
        {
            last = temp;
            last_best = true;
        }
        else
        {
            current = temp;
            last_best = false;
        }
    }
    return dir_best;
}

template <typename S>
static inline S _bisectionDistance(const MinkowskiDiff<S> & shape, Vector3<S> & dir, int max_iterations, S tol, Vector3<S> & p1, Vector3<S> & p2)
{ // bisection search method
    Vector3<S> last = shape.support(dir); // last support point
    S dist_best = dir.dot(last);
    for (unsigned int iterations = 0u; iterations < max_iterations; ++iterations)
    {
        Vector3<S> normal; // steepest rotation direction
        if (_bisectionOptimal(dir, last, last.norm() * tol, normal)) break;
      
        Vector3<S> current;
        bool last_best = true;
        bool optimal = false;
        for (; iterations < max_iterations; ++iterations) // find bisection bounds
        {
            dir = (-last).normalized();
            current = shape.support(dir);
            S dist = dir.dot(current);
            if (dist > dist_best) break;
            dist_best = dist;
            Vector3<S> tnormal;
            if (_bisectionOptimal(dir, current, current.norm() * tol, tnormal)) {
                optimal = true;
                break;
            }
            if (tnormal.dot(normal) < 0.0)
            {
                last_best = false;
                break;
            }
            last = current;
            normal = tnormal;
        }
        if (optimal) break;
        dir = _bisectionDistanceCore(shape, max_iterations, tol, iterations, dist_best, last_best, last, current, normal, p1, p2);
    }

    return dist_best;
}

// Computes the distance between two non-penetrating convex objects, returning the distance and
// nearest points on each object.
// @param p1 If the objects are non-penetrating, the point on the surface of
// obj1 closest to obj2 (expressed in the world frame).
// @param p2 If the objects are non-penetrating, the point on the surface of
// obj2 closest to obj1 (expressed in the world frame).
// @returns The minimum distance between the two objects. If they are
// penetrating, -1 is returned.
// @note Unlike _ccdDist function, this function does not need a warm-started
// simplex as the input argument.
template <typename S>
static inline S ccdGJKDist2(const MinkowskiDiff<S>& shape, unsigned int max_iterations, S tol, Vector3<S> &p1, Vector3<S> &p2)
{
  Simplex<S> simplex;
  // first find an intersection
  if (__ccdGJK(shape, simplex, max_iterations) == 0)
    return -1.0;
  return _ccdDist(shape, simplex, max_iterations, tol, p1, p2);
}

template <typename S>
static inline S ccdGJKDist2S(const MinkowskiDiff<S>& shape, unsigned int max_iterations, S tol, Vector3<S> &p1, Vector3<S> &p2)
{
  Simplex<S> simplex;
  // first find an intersection
  if (__ccdGJK(shape, simplex, max_iterations) == 0)
    return -1.0;
  S last_dist = std::numeric_limits<S>::max();
  Vector3<S> closest_p;
  _simplexClosestP(simplex, closest_p, last_dist);
  Vector3<S> dir = (-closest_p).normalized();
  return -_bisectionDistance(shape, dir, max_iterations, tol, p1, p2);
}

} // namespace libccd_extension

template <typename S>
bool GJKCollide(const MinkowskiDiff<S>& /*shape*/, unsigned int /*max_iterations*/, S /*tolerance*/,
                Vector3<S>* /*contact_points*/, S* /*penetration_depth*/, Vector3<S>* /*normal*/)
{
    return true;
}

template <typename S>
bool SeparationDistance(const MinkowskiDiff<S>& shape, unsigned int max_iterations, S tolerance,
                 S *res, Vector3<S> *p1, Vector3<S> *p2)
{
  Vector3<S> p1_ = Vector3<S>::Zero(), p2_ = Vector3<S>::Zero();
  // NOTE(JS): p1_ and p2_ are set to zeros in order to suppress uninitialized
  // warning. It seems the warnings occur since libccd_extension::ccdGJKDist2
  // conditionally set p1_ and p2_. If this wasn't intentional then please
  // remove the initialization of p1_ and p2_, and change the function
  // libccd_extension::ccdGJKDist2(...) to always set p1_ and p2_.
  S dist = libccd_extension::ccdGJKDist2S(shape, max_iterations, tolerance, p1_, p2_);
  if (p1) *p1 = p1_;
  if (p2) *p2 = p2_;
  if (res) *res = dist;
  if (dist < 0.0)
    return false;
  return true;
}

template <typename S>
bool SignedDistance(const MinkowskiDiff<S>& /*shape*/, unsigned int /*max_iterations*/, S /*tolerance*/,
        S* /*res*/, Vector3<S>* /*p1*/, Vector3<S>* /*p2*/)
{
    return true;
}
/*
{
  Vector3<S> p1_ = Vector3<S>::Zero(), p2_ = Vector3<S>::Zero();
  // NOTE(JS): p1_ and p2_ are set to zeros in order to suppress uninitialized
  // warning. It seems the warnings occur since libccd_extension::ccdGJKDist2
  // conditionally set p1_ and p2_. If this wasn't intentional then please
  // remove the initialization of p1_ and p2_, and change the function
  // libccd_extension::ccdGJKDist2(...) to always set p1_ and p2_.
  S dist = libccd_extension::ccdGJKSignedDist(shape, max_iterations, tolerance, p1_, p2_);
  if (p1) *p1 = p1_;
  if (p2) *p2 = p2_;
  if (res) *res = dist;
  if (dist < 0.0)
    return false;
  else
    return true;
}
*/
} // namespace detail
} // namespace fcl

#endif
