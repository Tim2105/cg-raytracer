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

/// \class Camera Camera.h
/// The camera class stores the following camera parameters in its
/// member variables: eye, center, up, fovy, width, and height. The function primary_ray(x,y)
/// produces a ray from the camera center through the pixel (x,y).
class Camera
{
public:
    /// Default (empty) constructor
    Camera() {}

    /// Construct camera by specifying all required data. Calls init().
    /// \param[in] eye Camera center
    /// \param[in] center Center of the scene, which the camera is looking at.
    /// \param[in] up Vector specifying the up direction.
    /// \param[in] fovy Opening angle (field of view) in y-direction
    /// \param[in] width Width of the image (in pixels)
    /// \param[in] height Height of the image (in pixels)
    Camera(const vec3& eye, const vec3& center, const vec3& up, double fovy,
           unsigned int width, unsigned int height)
        : eye_(eye),
          center_(center),
          up_(up),
          fovy_(fovy),
          width_(width),
          height_(height)
    {
        init();
    }

    /// This function precomputes some variables that are later required for primary_ray()
    void init()
    {
        // compute viewing direction and distance of eye to scene center
        vec3 view = normalize(center_ - eye_);
        double dist = distance(center_, eye_);

        // compute width & height of the image plane
        // based on the opening angle of the camera (fovy) and the distance
        // of the eye to the near plane (dist)
        double w = width_;
        double h = height_;
        double image_height = 2.0 * dist * tan(0.5 * fovy_ / 180.0 * M_PI);
        double image_width = w / h * image_height;

        // compute right and up vectors on the image plane
        x_dir_ = normalize(cross(view, up_)) * image_width / w;
        y_dir_ = normalize(cross(x_dir_, view)) * image_height / h;

        // compute lower left corner on the image plane
        lower_left_ = center_ - 0.5f * w * x_dir_ - 0.5f * h * y_dir_;
    }

    /// create a ray for a pixel in the image
    /// \param[in] x pixel location in image
    /// \param[in] y pixel location in image
    Ray primary_ray(double x, double y) const
    {
        return Ray(eye_, lower_left_ + x * x_dir_ + y * y_dir_ - eye_);
    }

public:
    /// position of the eye in 3D space (camera center)
    vec3 eye_;

    /// the center of the scene the camera is looking at
    vec3 center_;

    /// up-direction of the camera
    vec3 up_;

    /// opening angle (field of view) in y-direction
    double fovy_;

    /// image width in pixels
    int width_;

    /// image height in pixels
    int height_;

private:
    vec3 x_dir_;
    vec3 y_dir_;
    vec3 lower_left_;
};

//-----------------------------------------------------------------------------

/// read camera data from stream
inline std::istream& operator>>(std::istream& is, Camera& c)
{
    is >> c.eye_ >> c.center_ >> c.up_ >> c.fovy_ >> c.width_ >> c.height_;
    c.init();
    return is;
}

//=============================================================================
