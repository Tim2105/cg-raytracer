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

    vec3 normal = normalize(normal_);
    double normalDirDot = dot(normal, ray.direction_);

    if(normalDirDot == 0.0)
        return false;

    double distToOrigin = dot(normal, center_);
    intersection_distance = (-dot(normal, ray.origin_) + distToOrigin) /
                              normalDirDot;

    if(intersection_distance <= 1e-5)
        return false;

    intersection_normal = normal;
    intersection_point = ray(intersection_distance);

    return true;
}

//=============================================================================
