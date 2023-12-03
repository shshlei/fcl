// old
template <typename S>
static inline S _bisectionDistance(const MinkowskiDiff<S>& shape, Simplex<S> &simplex, unsigned int max_iterations, S tol,
        Vector3<S> &p1, Vector3<S> &p2)
{
    S best_dist = -std::numeric_limits<S>::max(), dist = best_dist, distp = -best_dist;
    Vector3<S> dir; // direction vector
    Support<S> last, best_support; // last support point
    
    unsigned int iterations = 0u;
    Vector3<S> closest_p; 
    while (iterations < max_iterations)
    {
        ++iterations; 
        _simplexClosestP(simplex, closest_p, distp);
        _ccdSupport(shape, closest_p, dir, last);
        if (_ccdOptimal(closest_p, dir, last, tol)) // termination one 
        {
            extractClosestPoints(simplex, p1, p2, closest_p);
            return distp;
        }
        best_dist = dist = -dir.dot(last.v);
        //if (dist > 0.0)
            break;
        simplex.add(last);
    }
    best_support = last;
    Vector3<S> normal = best_support.v.cross(dir), tnormal;
    S len = normal.norm();
    if (len < dist * tol)
    {
        extractClosestPoints(simplex, p1, p2, closest_p);
        return distp;
    }
    normal /= len;
    while (iterations < max_iterations)
    {
        ++iterations; 
        while (true)
        {
            _ccdSupport(shape, best_support.v, dir, last);
            dist = -dir.dot(last.v);
            if (dist <= best_dist)
                break;
            tnormal = last.v.cross(dir);
            S len = tnormal.norm();
            if (len < dist * tol) // termination
            {
                p1 = shape.support0(dir).dot(dir)*dir;
                p2 = shape.support1(-dir).dot(dir)*dir;
                return dist;
            }
            tnormal /= len;
            if (dist - best_dist < tol && tnormal.dot(normal) < 0.0)
                break;
            best_dist = dist;
            best_support = last;
            normal = tnormal;
        }
        Vector3<S> temp1 = best_support.v, temp2 = last.v, bnormal = normal;
        if (dist > best_dist)
        {
            best_dist = dist;
            best_support = last;
            bnormal = tnormal;
        }
        S last_dist = best_dist;
        while (true) // bisection core
        {
            dir = normal.cross(temp1 - temp2).normalized();
            last = shape.supportS(dir);
            dist = -dir.dot(last.v);
            tnormal = last.v.cross(dir);
            if (dist >= best_dist)
            {
                best_support = last;
                best_dist = dist;
                S len = tnormal.norm();
                if (len < dist * tol) // termination
                {
                    p1 = shape.support0(dir).dot(dir)*dir;
                    p2 = shape.support1(-dir).dot(dir)*dir;
                    return dist;
                }
                bnormal = tnormal / len;
            }
            if (_ccdOptimal(temp1, dir, last, tol))
                break;
            S d = tnormal.dot(normal);
            if (d >= 0.0)
                temp1 = last.v;
            else 
                temp2 = last.v;
        }
        if (best_dist - last_dist < tol)
            break;
        normal = bnormal;
    }

    /*
    simplex.set(0, best_support);
    simplex.setSize(1);
    while (iterations < max_iterations)
    {
        ++iterations; 
        _simplexClosestP(simplex, closest_p, distp);
        _ccdSupport(shape, closest_p, dir, last);
        if (_ccdOptimal(closest_p, dir, last, tol)) // termination one 
        {
            extractClosestPoints(simplex, p1, p2, closest_p);
            return distp;
        }
        simplex.add(last);
    }
    */

    return best_dist;
}
/*
{
    S best_dist = -std::numeric_limits<S>::max(), dist = best_dist, distp = -best_dist;
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
            if (_ccdOptimal(closest_p, dir, last, tol)) // termination one 
            {
                extractClosestPoints(simplex, p1, p2, closest_p);
                return distp;
            }
        }

        dist = -dir.dot(last.v);
        bool invalid = std::max(dist, best_dist) <= 0.0;
        if (invalid || iterations == 0u)
        {
            simplified = true;
            if (dist > best_dist)
            {
                best_dist = dist;
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
        if (dist > best_dist)
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
            if (simplified)
            {
                S last_dist = right ? dist : best_dist;
                Support<S> temp1, temp2;
                if (right)
                {
                    temp1 = last;
                    temp2 = best_support;
                }
                else 
                {
                    temp1 = best_support;
                    temp2 = last;
                }
                Vector3<S> normal = temp1.v.cross(p_B);
                S len = normal.norm();
                if (len < last_dist * tol) // termination
                {
                    p1 = shape.support0(p_B).dot(p_B)*p_B;
                    p2 = shape.support1(-p_B).dot(p_B)*p_B;
                    return last_dist;
                }
                normal /= len;
                if (right)
                    bisection = temp2.v.cross(best_dir).dot(normal) < 0.0;
                else 
                    bisection = temp2.v.cross(dir).dot(normal) < 0.0;
                if (bisection)
                {
                    while (true)
                    {
                        Vector3<S> dir1 = normal.cross(temp1.v - temp2.v).normalized();
                        Support<S> temp = shape.supportS(dir1);
                        S d = -dir1.dot(temp.v);
                        Vector3<S> tnormal = temp.v.cross(dir1);
                        if (d > last_dist)
                        {
                            dist = d;
                            dir = dir1; last = temp;
                            if (dist - last_dist < 10.0 * tol) // todo
                                break;
                            last_dist = dist;

                            S len = tnormal.norm();
                            if (len < d * tol) // termination
                            {
                                p1 = shape.support0(dir1).dot(dir1)*dir1;
                                p2 = shape.support1(-dir1).dot(dir1)*dir1;
                                return d;
                            }
                        }
                        if (_ccdOptimal(temp1.v, dir1, temp, tol))
                            break;
                        d = tnormal.dot(normal);
                        if (d >= 0.0)
                            temp1 = temp;
                        else 
                            temp2 = temp;
                    }
                }
            }
            else 
            {
                S last_dist = right ? dist : best_dist;
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
                        s *= 0.5;
                    }
                }
            }
        }

        if (dist >= best_dist && dist - best_dist < tol) // termination two // todo
        {
            p1 = shape.support0(dir).dot(dir)*dir;
            p2 = shape.support1(-dir).dot(dir)*dir;
            return dist;
        }
        if (dist > best_dist)
        {
            simplified = true;
            best_dir = dir;
            best_support = last;
            best_dist = dist;
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

    return best_dist;
}
*/

