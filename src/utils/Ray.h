//=============================================================================
//
//   Exercise code for the lecture
//   "Computer Graphics"
//   by Prof. Dr. Mario Botsch, TU Dortmund
//
//   Copyright (C) Computer Graphics Group, TU Dortmund.
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include "vec3.h"

//== CLASS DEFINITION =========================================================

/// \class Ray Ray.h
/// This class implements a ray, specified by its origin and direction.
/// It provides a convenient function to compute the point ray(t) at a specific
/// ray paramter t.
class Ray
{
public:
    /// Constructor with origin and direction. Direction will be normalized.
    Ray(const vec3& origin = vec3(0, 0, 0),
        const vec3& direction = vec3(0, 0, 1))
    {
        origin_ = origin;
        direction_ = normalize(direction); // normalize direction
    }

    /// Compute the point on the ray at the parameter \c t, which is
    /// origin + t*direction.
    vec3 operator()(double t) const { return origin_ + t * direction_; }

public:
    /// origin of the ray
    vec3 origin_;
    /// direction of the ray (should be normalized)
    vec3 direction_;
};

//-----------------------------------------------------------------------------

/// read ray from stream
inline std::istream& operator>>(std::istream& is, Ray& r)
{
    is >> r.origin_ >> r.direction_;
    return is;
}

//=============================================================================
