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

/** \author Jia Pan */

#include <gtest/gtest.h>

#include "fcl/traversal/traversal_node_bvhs.h"
#include "fcl/traversal/traversal_node_setup.h"
#include "fcl/collision_node.h"
#include "fcl/collision.h"
#include "fcl/BV/BV.h"
#include "fcl/shape/geometric_shapes.h"
#include "fcl/narrowphase/narrowphase.h"
#include "test_fcl_utility.h"
#include "fcl_resources/config.h"

using namespace fcl;

template<typename BV>
bool collide_Test(const Transform3d& tf,
                  const std::vector<Vector3d>& vertices1, const std::vector<Triangle>& triangles1,
                  const std::vector<Vector3d>& vertices2, const std::vector<Triangle>& triangles2, SplitMethodType split_method, bool verbose = true);

template<typename BV>
bool collide_Test2(const Transform3d& tf,
                   const std::vector<Vector3d>& vertices1, const std::vector<Triangle>& triangles1,
                   const std::vector<Vector3d>& vertices2, const std::vector<Triangle>& triangles2, SplitMethodType split_method, bool verbose = true);

template<typename BV, typename TraversalNode>
bool collide_Test_Oriented(const Transform3d& tf,
                           const std::vector<Vector3d>& vertices1, const std::vector<Triangle>& triangles1,
                           const std::vector<Vector3d>& vertices2, const std::vector<Triangle>& triangles2, SplitMethodType split_method, bool verbose = true);


template<typename BV>
bool test_collide_func(const Transform3d& tf,
                       const std::vector<Vector3d>& vertices1, const std::vector<Triangle>& triangles1,
                       const std::vector<Vector3d>& vertices2, const std::vector<Triangle>& triangles2, SplitMethodType split_method);

int num_max_contacts = std::numeric_limits<int>::max();
bool enable_contact = true;

std::vector<Contact> global_pairs;
std::vector<Contact> global_pairs_now;

GTEST_TEST(FCL_COLLISION, OBB_Box_test)
{
  FCL_REAL r_extents[] = {-1000, -1000, -1000, 1000, 1000, 1000};
  std::vector<Transform3d> rotate_transform;
  generateRandomTransforms(r_extents, rotate_transform, 1);
  
  AABB aabb1;
  aabb1.min_ = Vector3d(-600, -600, -600);
  aabb1.max_ = Vector3d(600, 600, 600);

  OBBd obb1;
  convertBV(aabb1, rotate_transform[0], obb1);
  Boxd box1;
  Transform3d box1_tf;
  constructBox(aabb1, rotate_transform[0], box1, box1_tf);

  FCL_REAL extents[] = {-1000, -1000, -1000, 1000, 1000, 1000};
  std::size_t n = 1000;

  std::vector<Transform3d> transforms;
  generateRandomTransforms(extents, transforms, n);

  for(std::size_t i = 0; i < transforms.size(); ++i)
  {
    AABB aabb;
    aabb.min_ = aabb1.min_ * 0.5;
    aabb.max_ = aabb1.max_ * 0.5;    

    OBBd obb2;
    convertBV(aabb, transforms[i], obb2);
    
    Boxd box2;
    Transform3d box2_tf;
    constructBox(aabb, transforms[i], box2, box2_tf);

    GJKSolver_libccd solver;

    bool overlap_obb = obb1.overlap(obb2);
    bool overlap_box = solver.shapeIntersect(box1, box1_tf, box2, box2_tf, NULL);
    
    EXPECT_TRUE(overlap_obb == overlap_box);
  }
}

GTEST_TEST(FCL_COLLISION, OBB_shape_test)
{
  FCL_REAL r_extents[] = {-1000, -1000, -1000, 1000, 1000, 1000};
  std::vector<Transform3d> rotate_transform;
  generateRandomTransforms(r_extents, rotate_transform, 1);
  
  AABB aabb1;
  aabb1.min_ = Vector3d(-600, -600, -600);
  aabb1.max_ = Vector3d(600, 600, 600);

  OBBd obb1;
  convertBV(aabb1, rotate_transform[0], obb1);
  Boxd box1;
  Transform3d box1_tf;
  constructBox(aabb1, rotate_transform[0], box1, box1_tf);

  FCL_REAL extents[] = {-1000, -1000, -1000, 1000, 1000, 1000};
  std::size_t n = 1000;

  std::vector<Transform3d> transforms;
  generateRandomTransforms(extents, transforms, n);

  for(std::size_t i = 0; i < transforms.size(); ++i)
  {
    FCL_REAL len = (aabb1.max_[0] - aabb1.min_[0]) * 0.5;
    OBBd obb2;
    GJKSolver_libccd solver;
 
    {  
      Sphered sphere(len);
      computeBV(sphere, transforms[i], obb2);
 
      bool overlap_obb = obb1.overlap(obb2);
      bool overlap_sphere = solver.shapeIntersect(box1, box1_tf, sphere, transforms[i], NULL);
      EXPECT_TRUE(overlap_obb >= overlap_sphere);
    }

    {
      Ellipsoidd ellipsoid(len, len, len);
      computeBV(ellipsoid, transforms[i], obb2);

      bool overlap_obb = obb1.overlap(obb2);
      bool overlap_ellipsoid = solver.shapeIntersect(box1, box1_tf, ellipsoid, transforms[i], NULL);
      EXPECT_TRUE(overlap_obb >= overlap_ellipsoid);
    }

    {
      Capsuled capsule(len, 2 * len);
      computeBV(capsule, transforms[i], obb2);
      
      bool overlap_obb = obb1.overlap(obb2);
      bool overlap_capsule = solver.shapeIntersect(box1, box1_tf, capsule, transforms[i], NULL);
      EXPECT_TRUE(overlap_obb >= overlap_capsule);
    }

    {
      Coned cone(len, 2 * len);
      computeBV(cone, transforms[i], obb2);
      
      bool overlap_obb = obb1.overlap(obb2);
      bool overlap_cone = solver.shapeIntersect(box1, box1_tf, cone, transforms[i], NULL);
      EXPECT_TRUE(overlap_obb >= overlap_cone);
    }

    {
      Cylinderd cylinder(len, 2 * len);
      computeBV(cylinder, transforms[i], obb2);
      
      bool overlap_obb = obb1.overlap(obb2);
      bool overlap_cylinder = solver.shapeIntersect(box1, box1_tf, cylinder, transforms[i], NULL);
      EXPECT_TRUE(overlap_obb >= overlap_cylinder);
    }
  }
}

GTEST_TEST(FCL_COLLISION, OBB_AABB_test)
{
  FCL_REAL extents[] = {-1000, -1000, -1000, 1000, 1000, 1000};
  std::size_t n = 1000;

  std::vector<Transform3d> transforms;
  generateRandomTransforms(extents, transforms, n);

  AABB aabb1;
  aabb1.min_ = Vector3d(-600, -600, -600);
  aabb1.max_ = Vector3d(600, 600, 600);
  
  OBBd obb1;
  convertBV(aabb1, Transform3d::Identity(), obb1);
  
  for(std::size_t i = 0; i < transforms.size(); ++i)
  {
    AABB aabb;
    aabb.min_ = aabb1.min_ * 0.5;
    aabb.max_ = aabb1.max_ * 0.5;    

    AABB aabb2 = translate(aabb, transforms[i].translation());
    
    OBBd obb2;
    convertBV(aabb, Transform3d(Eigen::Translation3d(transforms[i].translation())), obb2);

    bool overlap_aabb = aabb1.overlap(aabb2);
    bool overlap_obb = obb1.overlap(obb2);
    if(overlap_aabb != overlap_obb)
    {
      std::cout << aabb1.min_.transpose() << " " << aabb1.max_.transpose() << std::endl;
      std::cout << aabb2.min_.transpose() << " " << aabb2.max_.transpose() << std::endl;
      std::cout << obb1.To.transpose() << " " << obb1.extent.transpose() << " " << obb1.axis.col(0).transpose() << " " << obb1.axis.col(1).transpose() << " " << obb1.axis.col(2).transpose() << std::endl;
      std::cout << obb2.To.transpose() << " " << obb2.extent.transpose() << " " << obb2.axis.col(0).transpose() << " " << obb2.axis.col(1).transpose() << " " << obb2.axis.col(2).transpose() << std::endl;
    }

    EXPECT_TRUE(overlap_aabb == overlap_obb);
  }
  std::cout << std::endl;
}

GTEST_TEST(FCL_COLLISION, mesh_mesh)
{
  std::vector<Vector3d> p1, p2;
  std::vector<Triangle> t1, t2;
  
  loadOBJFile(TEST_RESOURCES_DIR"/env.obj", p1, t1);
  loadOBJFile(TEST_RESOURCES_DIR"/rob.obj", p2, t2);

  std::vector<Transform3d> transforms;
  FCL_REAL extents[] = {-3000, -3000, 0, 3000, 3000, 3000};
#ifdef FCL_BUILD_TYPE_DEBUG
  std::size_t n = 1;
#else
  std::size_t n = 10;
#endif
  bool verbose = false;

  generateRandomTransforms(extents, transforms, n);

  // collision
  for(std::size_t i = 0; i < transforms.size(); ++i)
  {
    global_pairs.clear();
    global_pairs_now.clear();

    collide_Test<OBBd>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);

    collide_Test<OBBd>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<OBBd>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<RSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<RSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<RSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<AABB>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<AABB>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<AABB>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<KDOP<24> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<KDOP<24> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<KDOP<24> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<KDOP<18> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<KDOP<18> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<KDOP<18> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<KDOP<16> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<KDOP<16> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<KDOP<16> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<OBBd>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<OBBd>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<OBBd>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<RSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<RSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<RSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<AABB>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<AABB>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<AABB>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<KDOP<24> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<KDOP<24> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<KDOP<24> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<KDOP<18> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<KDOP<18> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<KDOP<18> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<KDOP<16> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<KDOP<16> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<KDOP<16> >(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test_Oriented<OBBd, MeshCollisionTraversalNodeOBB>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);

    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test_Oriented<OBBd, MeshCollisionTraversalNodeOBB>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);

    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test_Oriented<OBBd, MeshCollisionTraversalNodeOBB>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test_Oriented<RSS, MeshCollisionTraversalNodeRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test_Oriented<RSS, MeshCollisionTraversalNodeRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test_Oriented<RSS, MeshCollisionTraversalNodeRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    test_collide_func<RSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    test_collide_func<OBBd>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    test_collide_func<AABB>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    
    collide_Test<kIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<kIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<kIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<kIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<kIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<kIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }
    
    collide_Test_Oriented<kIOS, MeshCollisionTraversalNodekIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test_Oriented<kIOS, MeshCollisionTraversalNodekIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test_Oriented<kIOS, MeshCollisionTraversalNodekIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    test_collide_func<kIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    test_collide_func<kIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    test_collide_func<kIOS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<OBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<OBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test<OBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<OBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<OBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test2<OBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }
    
    collide_Test_Oriented<OBBRSS, MeshCollisionTraversalNodeOBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test_Oriented<OBBRSS, MeshCollisionTraversalNodeOBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    collide_Test_Oriented<OBBRSS, MeshCollisionTraversalNodeOBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER, verbose);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    test_collide_func<OBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEAN);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    test_collide_func<OBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_MEDIAN);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }

    test_collide_func<OBBRSS>(transforms[i], p1, t1, p2, t2, SPLIT_METHOD_BV_CENTER);
    EXPECT_TRUE(global_pairs.size() == global_pairs_now.size());
    for(std::size_t j = 0; j < global_pairs.size(); ++j)
    {
      EXPECT_TRUE(global_pairs[j].b1 == global_pairs_now[j].b1);
      EXPECT_TRUE(global_pairs[j].b2 == global_pairs_now[j].b2);
    }
  }
}


template<typename BV>
bool collide_Test2(const Transform3d& tf,
                   const std::vector<Vector3d>& vertices1, const std::vector<Triangle>& triangles1,
                   const std::vector<Vector3d>& vertices2, const std::vector<Triangle>& triangles2, SplitMethodType split_method, bool verbose)
{
  BVHModel<BV> m1;
  BVHModel<BV> m2;
  m1.bv_splitter.reset(new BVSplitter<BV>(split_method));
  m2.bv_splitter.reset(new BVSplitter<BV>(split_method));

  std::vector<Vector3d> vertices1_new(vertices1.size());
  for(unsigned int i = 0; i < vertices1_new.size(); ++i)
  {
    vertices1_new[i] = tf * vertices1[i];
  }


  m1.beginModel();
  m1.addSubModel(vertices1_new, triangles1);
  m1.endModel();

  m2.beginModel();
  m2.addSubModel(vertices2, triangles2);
  m2.endModel();

  Transform3d pose1 = Transform3d::Identity();
  Transform3d pose2 = Transform3d::Identity();

  CollisionResult local_result;
  MeshCollisionTraversalNode<BV> node;

  if(!initialize<BV>(node, m1, pose1, m2, pose2,
                     CollisionRequest(num_max_contacts, enable_contact), local_result))
    std::cout << "initialize error" << std::endl;

  node.enable_statistics = verbose;

  collide(&node);


  if(local_result.numContacts() > 0)
  {
    if(global_pairs.size() == 0)
    {
      local_result.getContacts(global_pairs);
      std::sort(global_pairs.begin(), global_pairs.end());
    }
    else
    {
      local_result.getContacts(global_pairs_now);
      std::sort(global_pairs_now.begin(), global_pairs_now.end());
    }


    if(verbose)
      std::cout << "in collision " << local_result.numContacts() << ": " << std::endl;
    if(verbose) std::cout << node.num_bv_tests << " " << node.num_leaf_tests << std::endl;
    return true;
  }
  else
  {
    if(verbose) std::cout << "collision free " << std::endl;
    if(verbose) std::cout << node.num_bv_tests << " " << node.num_leaf_tests << std::endl;
    return false;
  }
}

template<typename BV>
bool collide_Test(const Transform3d& tf,
                  const std::vector<Vector3d>& vertices1, const std::vector<Triangle>& triangles1,
                  const std::vector<Vector3d>& vertices2, const std::vector<Triangle>& triangles2, SplitMethodType split_method, bool verbose)
{
  BVHModel<BV> m1;
  BVHModel<BV> m2;
  m1.bv_splitter.reset(new BVSplitter<BV>(split_method));
  m2.bv_splitter.reset(new BVSplitter<BV>(split_method));

  m1.beginModel();
  m1.addSubModel(vertices1, triangles1);
  m1.endModel();

  m2.beginModel();
  m2.addSubModel(vertices2, triangles2);
  m2.endModel();

  Transform3d pose1(tf);
  Transform3d pose2 = Transform3d::Identity();

  CollisionResult local_result;
  MeshCollisionTraversalNode<BV> node;

  if(!initialize<BV>(node, m1, pose1, m2, pose2,
                     CollisionRequest(num_max_contacts, enable_contact), local_result))
    std::cout << "initialize error" << std::endl;

  node.enable_statistics = verbose;

  collide(&node);


  if(local_result.numContacts() > 0)
  {
    if(global_pairs.size() == 0)
    {
      local_result.getContacts(global_pairs);
      std::sort(global_pairs.begin(), global_pairs.end());
    }
    else
    {
      local_result.getContacts(global_pairs_now);
      std::sort(global_pairs_now.begin(), global_pairs_now.end());
    }

    if(verbose)
      std::cout << "in collision " << local_result.numContacts() << ": " << std::endl;
    if(verbose) std::cout << node.num_bv_tests << " " << node.num_leaf_tests << std::endl;
    return true;
  }
  else
  {
    if(verbose) std::cout << "collision free " << std::endl;
    if(verbose) std::cout << node.num_bv_tests << " " << node.num_leaf_tests << std::endl;
    return false;
  }
}

template<typename BV, typename TraversalNode>
bool collide_Test_Oriented(const Transform3d& tf,
                           const std::vector<Vector3d>& vertices1, const std::vector<Triangle>& triangles1,
                           const std::vector<Vector3d>& vertices2, const std::vector<Triangle>& triangles2, SplitMethodType split_method, bool verbose)
{
  BVHModel<BV> m1;
  BVHModel<BV> m2;
  m1.bv_splitter.reset(new BVSplitter<BV>(split_method));
  m2.bv_splitter.reset(new BVSplitter<BV>(split_method));

  m1.beginModel();
  m1.addSubModel(vertices1, triangles1);
  m1.endModel();

  m2.beginModel();
  m2.addSubModel(vertices2, triangles2);
  m2.endModel();

  Transform3d pose1(tf);
  Transform3d pose2 = Transform3d::Identity();

  CollisionResult local_result;
  TraversalNode node;
  if(!initialize(node, (const BVHModel<BV>&)m1, pose1, (const BVHModel<BV>&)m2, pose2, 
                 CollisionRequest(num_max_contacts, enable_contact), local_result))
    std::cout << "initialize error" << std::endl;

  node.enable_statistics = verbose;

  collide(&node);

  if(local_result.numContacts() > 0)
  {
    if(global_pairs.size() == 0)
    {
      local_result.getContacts(global_pairs);
      std::sort(global_pairs.begin(), global_pairs.end());
    }
    else
    {
      local_result.getContacts(global_pairs_now);
      std::sort(global_pairs_now.begin(), global_pairs_now.end());
    }

    if(verbose)
      std::cout << "in collision " << local_result.numContacts() << ": " << std::endl;
    if(verbose) std::cout << node.num_bv_tests << " " << node.num_leaf_tests << std::endl;
    return true;
  }
  else
  {
    if(verbose) std::cout << "collision free " << std::endl;
    if(verbose) std::cout << node.num_bv_tests << " " << node.num_leaf_tests << std::endl;
    return false;
  }
}


template<typename BV>
bool test_collide_func(const Transform3d& tf,
                       const std::vector<Vector3d>& vertices1, const std::vector<Triangle>& triangles1,
                       const std::vector<Vector3d>& vertices2, const std::vector<Triangle>& triangles2, SplitMethodType split_method)
{
  BVHModel<BV> m1;
  BVHModel<BV> m2;
  m1.bv_splitter.reset(new BVSplitter<BV>(split_method));
  m2.bv_splitter.reset(new BVSplitter<BV>(split_method));

  m1.beginModel();
  m1.addSubModel(vertices1, triangles1);
  m1.endModel();

  m2.beginModel();
  m2.addSubModel(vertices2, triangles2);
  m2.endModel();

  Transform3d pose1(tf);
  Transform3d pose2 = Transform3d::Identity();

  std::vector<Contact> contacts;

  CollisionRequest request(num_max_contacts, enable_contact);
  CollisionResult result;
  int num_contacts = collide(&m1, pose1, &m2, pose2, 
                             request, result);
	
  result.getContacts(contacts);

  global_pairs_now.resize(num_contacts);

  for(int i = 0; i < num_contacts; ++i)
  {
    global_pairs_now[i].b1 = contacts[i].b1;
    global_pairs_now[i].b2 = contacts[i].b2;
  }

  std::sort(global_pairs_now.begin(), global_pairs_now.end());

  if(num_contacts > 0) return true;
  else return false;
}

//==============================================================================
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
