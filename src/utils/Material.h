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
#include "Image.h"

//== CLASS DEFINITION =========================================================

/// \class Material Material.h
/// This class stores all material parameters.
struct Material
{
    /// ambient color
    vec3 ambient;

    /// diffuse color
    vec3 diffuse;

    /// specular color
    vec3 specular;

    /// shininess factor
    double shininess;

    /// reflectivity factor (1=perfect mirror, 0=no reflection).
    double mirror;

    /// determines if the material can be in shadow (important for sky material in pokemon scene)
    bool shadowable;

    /// Read in Texture png
    Image texture_png;
};

//-----------------------------------------------------------------------------

/// read material from stream
inline std::istream& operator>>(std::istream& is, Material& m)
{
    is >> m.ambient >> m.diffuse >> m.specular >> m.shininess >> m.mirror;
    if (m.diffuse[0] < 1e-5 && m.diffuse[1] < 1e-5 && m.diffuse[2] < 1e-5 &&
        m.specular[0] < 1e-5 && m.specular[1] < 1e-5 && m.specular[2] < 1e-5)
        m.shadowable = false;
    else
        m.shadowable = true;
    return is;
}

//=============================================================================
