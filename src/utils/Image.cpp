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

#include "Image.h"
#include <fstream>
// #include <lodepng.h>
#include <stb_image_write.h>
#include <stb_image.h>

//== CLASS DEFINITION =========================================================

Image::Image(unsigned int _width, unsigned int _height)
{
    resize(_width, _height);
}

void Image::resize(unsigned int _width, unsigned int _height)
{
    width_ = _width;
    height_ = _height;
    pixels_.resize(width_ * height_);
}

bool Image::write_png(const std::string& _filename)
{
    // convert float pixel in [0,1] values to unsigned char in [0,255]
    std::vector<unsigned char> rgb(3 * width_ * height_);
    for (unsigned int y = 0; y < height_; ++y)
    {
        for (unsigned int x = 0; x < width_; ++x)
        {
            for (unsigned int c = 0; c < 3; ++c)
            {
                rgb[3 * (y * width_ + x) + c] =
                    255.0 * pixels_[y * width_ + x][c];
            }
        }
    }

    // write RGB array to file
    stbi_flip_vertically_on_write(true);
    bool ok = stbi_write_png(_filename.c_str(), width_, height_, 3, rgb.data(),
                             3 * width_);
    return ok;
}

bool Image::read(const std::string& _filename)
{
    // load with stb_image
    int width, height, n;
    stbi_set_flip_vertically_on_load(true);
    unsigned char* rgb = stbi_load(_filename.c_str(), &width, &height, &n, 3);
    if (!rgb)
        return false;

    // allocate memory
    resize(width, height);

    // convert RGB (unsigned char) to float colors
    vec3 color;
    for (unsigned int y = 0; y < height_; ++y)
    {
        for (unsigned int x = 0; x < width_; ++x)
        {
            for (unsigned int c = 0; c < 3; ++c)
            {
                color[c] = rgb[3 * (y * width_ + x) + c] / 255.0;
            }
            (*this)(x, y) = color;
        }
    }

    // free memory
    stbi_image_free(rgb);

    return true;
}

//=============================================================================
