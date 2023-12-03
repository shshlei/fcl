#include "fcl/narrowphase/detail/bisection_solver.h"
#include "fcl/narrowphase/detail/gjk_solver_libccd.h"
#include "fcl/narrowphase/detail/gjk_solver_indep.h"
#include "test_fcl_utility.h"
#include <fstream>
#include <string>
#include <vector>

using namespace fcl;
using namespace detail;
using namespace test;

template<typename Shape1, typename Shape2>
void timeCollect(const Shape1 &s1, const aligned_vector<Transform3d> &transforms1,
        const Shape2 &s2, const aligned_vector<Transform3d> &transforms2,
        double &time1, double &time2, double &error1, double &time3, double &error2);

int main()
{
    // Box-Ellipsoid-Sphere-Cylinder-Capsule-Cone
    unsigned int num = 6, i = 0, j = 0;
    std::vector<std::vector<double>> vtime1(num, std::vector<double>(num, 0.0));
    std::vector<std::vector<double>> vtime2(num, std::vector<double>(num, 0.0));
    std::vector<std::vector<double>> vtime3(num, std::vector<double>(num, 0.0));
    std::vector<std::vector<double>> verror1(num, std::vector<double>(num, 0.0));
    std::vector<std::vector<double>> verror2(num, std::vector<double>(num, 0.0));

    for (unsigned int len = 0; len < 1; len++)
    {
        //srand((unsigned)time(NULL));  // todo

        aligned_vector<Transform3d> transforms1, transforms2;
        double extents[] = {-rand_interval(5.0, 100.0), -rand_interval(5.0, 100.0), -rand_interval(5.0, 100.0), rand_interval(5.0, 100.0), rand_interval(5.0, 100.0), rand_interval(5.0, 100.0)};
        std::size_t n = 1000;
        generateRandomTransforms(extents, transforms1, n);
        generateRandomTransforms(extents, transforms2, n);

        double low = rand_interval(0.1, 10.0), high = rand_interval(0.1, 10.0) + low;

        //==============================================================================
        //Box
        i = 0;
        {
            j = 0;
            Boxd s1(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            Boxd s2(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Boxd, Boxd>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Box-Box " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 << std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 1;
            Boxd s1(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            Ellipsoidd s2(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Boxd, Ellipsoidd>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Box-Ellipsoid " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 2;
            Boxd s1(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            Sphered s2(rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Boxd, Sphered>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Box-Sphere " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 3;
            Boxd s1(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            Cylinderd s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Boxd, Cylinderd>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Box-Cylinder " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 4;
            Boxd s1(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            Capsuled s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Boxd, Capsuled>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Box-Capsule " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 5;
            Boxd  s1(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            Coned s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Boxd, Coned>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Box-Cone " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        std::cout << std::endl << std::endl;

        //==============================================================================
        //Ellipsoid
        i = 1;
        {
            j = 1;
            Ellipsoidd s1(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            Ellipsoidd s2(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Ellipsoidd, Ellipsoidd>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Ellipsoid-Ellipsoidd " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 2;
            Ellipsoidd s1(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            Sphered s2(rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Ellipsoidd, Sphered>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Ellipsoid-Sphere " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 3;
            Ellipsoidd s1(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            Cylinderd s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Ellipsoidd, Cylinderd>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Ellipsoid-Cylinder " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 4;
            Ellipsoidd s1(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            Capsuled s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Ellipsoidd, Capsuled>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Ellipsoid-Capsule " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 5;
            Ellipsoidd s1(rand_interval(low, high), rand_interval(low, high), rand_interval(low, high));
            Coned s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Ellipsoidd, Coned>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Ellipsoid-Cone " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        std::cout << std::endl << std::endl;

        //==============================================================================
        //Sphere
        i = 2;
        {
            j = 2;
            Sphered s1(rand_interval(low, high));
            Sphered s2(rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Sphered, Sphered>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Sphere-Sphere " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 3;
            Sphered s1(rand_interval(low, high));
            Cylinderd s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Sphered, Cylinderd>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Sphere-Cylinder " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 4;
            Sphered s1(rand_interval(low, high));
            Capsuled s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Sphered, Capsuled>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Sphere-Capsule " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 5;
            Sphered s1(rand_interval(low, high));
            Coned s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Sphered, Coned>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Sphere-Cone " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }
        std::cout << std::endl << std::endl;

        //==============================================================================
        //Cylinder
        i = 3;
        {
            j = 3;
            Cylinderd s1(rand_interval(low, high), rand_interval(low, high));
            Cylinderd s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Cylinderd, Cylinderd>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Cylinder-Cylinder " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 4;
            Cylinderd s1(rand_interval(low, high), rand_interval(low, high));
            Capsuled s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Cylinderd, Capsuled>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Cylinder-Capsule " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 5;
            Cylinderd s1(rand_interval(low, high), rand_interval(low, high));
            Coned s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Cylinderd, Coned>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Cylinder-Cone " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }
        std::cout << std::endl << std::endl;

        //==============================================================================
        //Capsule
        i = 4;
        {
            j = 4;
            Capsuled s1(rand_interval(low, high), rand_interval(low, high));
            Capsuled s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Capsuled, Capsuled>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Capsule-Capsule " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }

        {
            j = 5;
            Capsuled s1(rand_interval(low, high), rand_interval(low, high));
            Coned s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Capsuled, Coned>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Capsule-Cone " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }
        std::cout << std::endl << std::endl;

        //==============================================================================
        //Cone
        i = 5;
        {
            j = 5;
            Coned s1(rand_interval(low, high), rand_interval(low, high));
            Coned s2(rand_interval(low, high), rand_interval(low, high));
            double time1 = 0.0, time2 = 0.0, error = 0.0, time3 = 0.0, error2 = 0.0;
            timeCollect<Coned, Coned>(s1, transforms1, s2, transforms2, time1, time2, error, time3, error2);
            std::cout << "Cone-Cone " << "libccd_time " << time1 << " indep_time " << time2 << " error " << error << " libccds_time " << time3 << " error2 " << error2 <<  std::endl;
            vtime1[i][j] += time1;
            vtime2[i][j] += time2;
            vtime3[i][j] += time3;
            verror1[i][j] += error;
            verror2[i][j] += error2;
        }
    }

    /*
    std::ofstream ofs1("vtime1.txt", std::ios::binary | std::ios::out);
    std::ofstream ofs2("vtime2.txt", std::ios::binary | std::ios::out);
    std::ofstream ofs3("vtime3.txt", std::ios::binary | std::ios::out);
    std::ofstream ofs4("verror1.txt", std::ios::binary | std::ios::out);
    std::ofstream ofs5("verror2.txt", std::ios::binary | std::ios::out);
    for (i = 0; i < num; i++)
    {
        for (j = 0; j < num; j++)
        {
            ofs1 << vtime1[i][j] << " ";
            ofs2 << vtime2[i][j] << " ";           
            ofs3 << vtime3[i][j] << " ";
            ofs4 << verror1[i][j] << " ";
            ofs5 << verror2[i][j] << " ";
        }
        ofs1 << std::endl;
        ofs2 << std::endl;
        ofs3 << std::endl;
        ofs4 << std::endl;
        ofs5 << std::endl;
    }
    ofs1.close();
    ofs2.close();
    ofs3.close();
    ofs4.close();
    ofs5.close();
    */

    return 0;
}

template<typename Shape1, typename Shape2>
void timeCollect(const Shape1 &s1, const aligned_vector<Transform3d> &transforms1,
        const Shape2 &s2, const aligned_vector<Transform3d> &transforms2,
        double &time1, double &time2, double &error1, double &time3, double &error2)
{
    time1 = time2 = time3 = error1 = error2 = 0.0;
    GJKSolver_libccdd gjk_solver_libccd;
    GJKSolver_indepd gjk_solver_indep;
    // BisectionSolverd bisection_solver;

    test::Timer timer_dis;
    for (std::size_t i = 0; i < transforms1.size(); i++)
    {
        for (std::size_t j = 0; j < transforms2.size(); j++)
        {
            double dist1, dist2, dist3;
            Vector3d p1, p2;

            timer_dis.start();
            bool res = gjk_solver_libccd.shapeDistance<Shape1, Shape2>(s1, transforms1[i], s2, transforms2[j], &dist1, &p1,  &p2);
            //bool res = bisection_solver.shapeDistanceCCD<Shape1, Shape2>(s1, transforms1[i], s2, transforms2[j], &dist1, &p1,  &p2);
            timer_dis.stop();

            if (res)
            {
                time1 += timer_dis.getElapsedTimeInSec();

                timer_dis.start();
                gjk_solver_indep.shapeDistance<Shape1, Shape2>(s1, transforms1[i], s2, transforms2[j], &dist2, &p1,  &p2);
                timer_dis.stop();
                time2 += timer_dis.getElapsedTimeInSec();
                error1 += abs(dist1 - dist2);

                timer_dis.start();
                // bisection_solver.shapeDistance<Shape1, Shape2>(s1, transforms1[i], s2, transforms2[j], &dist3, &p1,  &p2);
                timer_dis.stop();
                time3 += timer_dis.getElapsedTimeInSec();
                error2 += abs(dist1 - dist3);
            }
        }
    }
}
