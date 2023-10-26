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

#include "Plane.h"
#include <float.h>

//== CLASS DEFINITION =========================================================

Plane::Plane(const vec3& c, const vec3& n) : center_(c), normal_(n) {}

//-----------------------------------------------------------------------------

bool Plane::intersect(const Ray& ray, vec3& intersection_point,
                      vec3& intersection_normal, vec3& intersection_diffuse,
                      double& intersection_distance) const
{
    intersection_diffuse = material_.diffuse;

    /** \todo
 * - compute the intersection of the plane with `ray`
 * - if ray and plane are parallel there is no intersection
 * - otherwise compute intersection data and store it in `intersection_point`, `intersection_normal`, and `intersection_distance`.
 * - return whether there is an intersection for t>1e-5 (avoids "shadow acne").
*/


    return false;
}

//=============================================================================
