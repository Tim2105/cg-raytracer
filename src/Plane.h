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

/// \class Plane Plane.h
/// This class implements a simple plane object.
/// A plane is specified by a center point and a normal vector.
/// This class overrides the intersection method Object::intersect().
class Plane : public Object
{
public:
    /// constructor
    Plane(const vec3& center = vec3(0, 0, 0),
          const vec3& normal = vec3(0, 1, 0));

    /// Compute the intersection of the plane with \c _ray. Return whether
    /// there is an intersection. If there is one, return the intersection data.
    /// This function overrides Object::intersect().
    /// \param[in] ray the ray to intersect the plane with
    /// \param[out] intersection_point position of the intersection
    /// \param[out] intersection_normal normal vector at the intersection point
    /// \param[out] intersection_diffuse the diffuse color at intersection point
    /// \param[out] intersection_distance ray parameter at the intesection point
    virtual bool intersect(const Ray& _ray, vec3& intersection_point,
                           vec3& intersection_normal,
                           vec3& intersection_diffuse,
                           double& intersection_distance) const;

public:
    /// one (arbitrary) point on the plane
    vec3 center_;
    /// normal vector of the plane
    vec3 normal_;
    /// color of the plane
    vec3 color_;
};

//-----------------------------------------------------------------------------

/// read plane from stream
inline std::istream& operator>>(std::istream& is, Plane& p)
{
    is >> p.center_ >> p.normal_ >> p.material_;
    return is;
}

//=============================================================================
