//=============================================================================
//
//   Exercise code for the lecture
//   "Computer Graphics"
//   by Prof. Dr. Mario Botsch, TU Dortmund
//
//   Copyright (C) Computer Graphics Group, TU Dortmund.
//
//=============================================================================

//== INCLUDES =================================================================

#include "Sphere.h"
#include <float.h>

//== CLASS DEFINITION =========================================================

Sphere::Sphere(const vec3& c, double r) : center_(c), radius_(r) {}

//-----------------------------------------------------------------------------

bool Sphere::intersect(const Ray& ray, vec3& intersection_point,
                       vec3& intersection_normal, vec3& intersection_diffuse,
                       double& intersection_distance) const
{
    const vec3 dir = ray.direction_;
    const vec3 oc = ray.origin_ - center_;
    intersection_diffuse = material_.diffuse;

    double a = dot(dir, dir);
    double b = 2.0 * dot(dir, oc);
    double c = dot(oc, oc) - radius_ * radius_;
    double d = b * b - 4.0 * a * c;

    // is there an intersection?
    if (d >= 0.0)
    {
        d = sqrt(d);

        // the two intersection parameters
        double t1 = (-b - d) / (2.0 * a);
        double t2 = (-b + d) / (2.0 * a);

        // find closer intersection
        intersection_distance = DBL_MAX;
        if (t1 > 1e-5 && t1 < intersection_distance)
            intersection_distance = t1;
        if (t2 > 1e-5 && t2 < intersection_distance)
            intersection_distance = t2;

        // was the intersection not just a numerical problem?
        if (intersection_distance != DBL_MAX)
        {
            // return intersection data
            intersection_point = ray(intersection_distance);
            intersection_normal = (intersection_point - center_) / radius_;

            return true;
        }
    }

    return false;
}

//=============================================================================
