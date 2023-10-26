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

//== includes =================================================================

#define _MINWINDEF_ // Ansonsten wird auf Win10 irgendeine 'minwindef.h' inkludiert, die ein Makro 'min' definiert

#include "utils/Object.h"
#include "utils/Light.h"
#include "utils/Ray.h"
#include "utils/Material.h"
#include "utils/Image.h"
#include "utils/Camera.h"

//== CLASS DEFINITION =========================================================

/// \class Raytracer Raytracer.h
/// In this class the main raytracing routines are implemented.
class Raytracer
{
public:
    Raytracer(const std::string& filename);
    ~Raytracer();

    /// reads the scene to raytrace from the given `filename`
    void read_scene(const std::string& filename);

    /// main function where the image is allocated and the raytracing is started
    void compute_image();

    /// write png file of the resulting image
    void write_image(const std::string& filename);

private:
    /// Finds the closest colision point for the passed Ray and returns the computed color for this object.
    /**
    *	@param ray passed Ray
    *	@param depth holds the information, how many times the `ray` had been reflected. Goes from 0 to `max_depth_`. Should be used for recursive function call.
    **/
    vec3 trace(const Ray& ray, int depth);

    /// Computes the closest intersection point between a ray and all objects in the scene.
    /**
    *   	@param ray Ray that should be tested for intersections with all objects in the scene.
    *   	@param intersection_material returns material of intersected object.
    *   	@param intersection_point returns intersection point
    *   	@param intersection_normal returns normal at `_point`
    *   	@param intersection_distance returns distance between the `ray`'s origin and `point`
    *   	@return returns `true`, if there is an intersection point between `ray` and at least one object in the scene.
    **/
    bool intersect_scene(const Ray& ray, Material& intersection_material,
                         vec3& intersection_point, vec3& intersection_normal,
                         double& intersection_distance);

    /**
     * @brief Überprüft, ob der Weg zwischen einem Punkt und einer Lichtquelle
     * blockiert weg, d.h. ob der Punkt im Schatten der Lichtquelle liegt.
     * 
     * @param point Der, zu überprüfende, Punkt.
     * @param light Die Positon der Lichtquelle.
     * @return Ein Wahrheitswert der angibt, ob der Punkt im Schatten liegt.
     */
    bool shadow(const vec3& point, const vec3& light);

    /// Computes the phong lighting.
    /**
    *	@param point the point, whose color should be determined.
    *	@param normal `point`'s normal
    *	@param view normalized direction from the `point' to the viewer's position.
    * 	@param material holds material parameters of the `point`, that should be lit.
    */
    vec3 lighting(const vec3& point, const vec3& normal, const vec3& view,
                  const Material& material);

    //-----------------------------------------------------------------------------
    // Member Variables
    //-----------------------------------------------------------------------------

    /// camera stores eye position, view direction, and can generate primary rays
    Camera camera_;

    /// our image is just a large array of color values
    Image image_;

    /// array for all lights in the scene
    std::vector<Light> lights_;

    /// array for all the objects in the scene
    std::vector<Object_ptr> objects_;

    /// max recursion depth for mirroring
    int max_depth_;

    /// background color
    vec3 background_;

    /// global ambient light
    vec3 ambience_;
};

//=============================================================================
