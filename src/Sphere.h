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

#include "utils/Object.h"
#include "utils/vec3.h"

//== CLASS DEFINITION =========================================================

/// \class Sphere Sphere.h
/// This class implements a sphere object, which is specified by its center
/// and its radius. This class overrides the intersection method Object::intersect().
class Sphere : public Object
{
public:
    /// Construct a sphere by specifying center and radius
    Sphere(const vec3& center = vec3(0, 0, 0), double radius = 1);

    /// Compute the intersection of the sphere with \c _ray. Return whether
    /// there is an intersection. If there is one, return the intersection data.
    /// This function overrides Object::intersect().
    /// \param[in] ray the ray to intersect the plane with
    /// \param[out] intersection_point position of the intersection
    /// \param[out] intersection_normal normal vector at the intersection point
    /// \param[out] intersection_diffuse the diffuse color at intersection point
    /// \param[out] intersection_distance ray parameter at the intesection point
    virtual bool intersect(const Ray& ray, vec3& intersection_point,
                           vec3& intersection_normal,
                           vec3& intersection_diffuse,
                           double& intersection_distance) const;

public:
    /// center position of the sphere
    vec3 center_;

    /// radius of the sphere
    double radius_;
};

//-----------------------------------------------------------------------------

/// read sphere from stream
inline std::istream& operator>>(std::istream& is, Sphere& s)
{
    is >> s.center_ >> s.radius_ >> s.material_;
    return is;
}

//=============================================================================
