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

#include "Ray.h"
#include "vec3.h"
#include "Material.h"

//== CLASS DEFINITION =========================================================

/// \class Object Object.h
/// This class implements an abstract class for an object.
/// Every derived object type will inherit the material property, and it
/// will have to override the virtual function Object::intersect().
struct Object
{
public:
    /// Default (empty) constructor
    Object() {}

    /// Destructor (if a class has virtual functions, then its destructor
    /// has to be virtual as well).
    virtual ~Object() {}

    /// Intersect the object with \c _ray, return whether there is an intersection.
    /// If \c _ray intersects the object, provide the following results:
    /// \param[in] ray the ray to intersect the object with
    /// \param[out] intersection_point the point of intersection
    /// \param[out] intersection_normal the surface normal at intersection point
    /// \param[out] intersection_diffuse the diffuse color at intersection point
    /// \param[out] intersection_t ray parameter at intersection point
    virtual bool intersect(const Ray& ray, vec3& intersection_point,
                           vec3& intersection_normal,
                           vec3& intersection_diffuse,
                           double& intersection_distance) const = 0;

    /// The material of this object
    Material material_;
};

/// This typedef is for convencience only
typedef Object* Object_ptr;

//=============================================================================
