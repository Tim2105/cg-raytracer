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

/// \class Light Light.h
/// The class represents a light source and stores position and color.
struct Light
{
    /// position of the light source
    vec3 position;

    /// color of the light source
    vec3 color;
};

//-----------------------------------------------------------------------------

/// read camera data from stream
inline std::istream& operator>>(std::istream& is, Light& l)
{
    is >> l.position >> l.color;
    return is;
}

//=============================================================================
