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
#include <vector>
#include <assert.h>

//== CLASS DEFINITION =========================================================

/// \class Image Image.h
/// This class stores an image as a big array of colors, which are represented
/// as vec3 values for RGB.
class Image
{
public:
    /// Construct an image of size _width times _height
    /// \param _width Width of the image in pixels
    /// \param _height Height of the image in pixels
    Image(unsigned int _width = 0, unsigned int _height = 0);

    /// Resize the image
    /// \param _width New width of the image in pixels
    /// \param _height New height of the image in pixels
    void resize(unsigned int _width, unsigned int _height);

    /// Returns image width in pixels.
    unsigned int width() const { return width_; }

    /// Returns image height in pixels.
    unsigned int height() const { return height_; }

    /// Read/write access to pixel (_x,_y). Use this to set the color by
    /// image(x,y) = color;
    vec3& operator()(unsigned int _x, unsigned int _y)
    {
        assert(_x < width_);
        assert(_y < height_);
        unsigned int el = _y * static_cast<unsigned int>(width_) + _x;
        return pixels_[el];
    }

    /// Read access to pixel (_x,_y).
    const vec3& operator()(unsigned int _x, unsigned int _y) const
    {
        assert(_x < width_);
        assert(_y < height_);
        return pixels_[_y * static_cast<unsigned int>(width_) + _x];
    }

    /// Writes the image in PNG format to a file.
    /// \param[in] _filename Filename to save the image to.
    bool write_png(const std::string& _filename);

    /// Reads an image out of a file in PNG format.
    /// \param[in] _filename Filename to read the image off.
    bool read(const std::string& _filename);

private:
    /// vector with all pixels in the image
    std::vector<vec3> pixels_;

    /// image width in pixels
    unsigned short int width_;

    /// image height in pixels
    unsigned short int height_;
};

//=============================================================================
