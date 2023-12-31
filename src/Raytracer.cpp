//=============================================================================
//
//   Exercise code for the lecture
//   "Computer Graphics"
//   by Prof. Dr. Mario Botsch, TU Dortmund
//
//   Copyright (C) Computer Graphics Group, TU Dortmund.
//
//=============================================================================

//== includes =================================================================

#include "Raytracer.h"

#include "utils/StopWatch.h"
#include "Sphere.h"
#include "Plane.h"
#include "Mesh.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <float.h>

//-----------------------------------------------------------------------------

Raytracer::Raytracer(const std::string& filename)
{
    read_scene(filename);
}

//-----------------------------------------------------------------------------

Raytracer::~Raytracer()
{
    //clean up
    for (auto o : objects_)
    {
        delete o;
    }
}

//-----------------------------------------------------------------------------

void Raytracer::read_scene(const std::string& filename)
{
    //clean up
    for (auto o : objects_)
    {
        delete o;
    }
    objects_.clear();
    lights_.clear();

    std::ifstream ifs(filename);
    if (!ifs)
    {
        std::cerr << "Cannot open file " << filename << std::endl;
        exit(1);
    }

    char line[200];
    std::string token;

    // parse file
    while (ifs && (ifs >> token) && (!ifs.eof()))
    {
        if (token[0] == '#')
        {
            ifs.getline(line, 200);
        }
        else if (token == "depth")
        {
            ifs >> max_depth_;
        }
        else if (token == "camera")
        {
            ifs >> camera_;
        }
        else if (token == "background")
        {
            ifs >> background_;
        }
        else if (token == "ambience")
        {
            ifs >> ambience_;
        }
        else if (token == "light")
        {
            Light light;
            ifs >> light;
            lights_.push_back(light);
        }
        else if (token == "plane")
        {
            Plane* p = new Plane;
            ifs >> (*p);
            objects_.push_back(p);
        }
        else if (token == "sphere")
        {
            Sphere* sphere = new Sphere();
            ifs >> (*sphere);
            objects_.push_back(sphere);
        }
        else if (token == "mesh")
        {
            std::string fn, mode;
            ifs >> fn >> mode;

            // add path of scene-file to mesh's filename
            std::string path(filename);
            path = path.substr(0, path.find_last_of('/') + 1);
            fn = path + fn;

            Mesh* mesh =
                new Mesh((mode == "FLAT" ? Mesh::FLAT : Mesh::PHONG), fn);

            ifs >> mesh->material_;

            objects_.push_back(mesh);
        }
    }
    ifs.close();

    std::cout << "\ndone (" << objects_.size() << " objects)\n";
}

//-----------------------------------------------------------------------------

#include <functional>

template <size_t n>
void genOffsets(double pixelOffsets[n][n][2]) {
    double pixelStepSize = 1.0 / n;

    double xOffset = -0.5 + pixelStepSize / 2;
    for(size_t i = 0; i < n; i++) {
        double yOffset = -0.5 + pixelStepSize / 2;
        for(size_t j = 0; j < n; j++) {
            pixelOffsets[i][j][0] = xOffset;
            pixelOffsets[i][j][1] = yOffset;

            yOffset += pixelStepSize;
        }
        xOffset += pixelStepSize;
    }
}

// template <size_t n>
// void ssTrace()

void Raytracer::compute_image()
{
    // allocate memory by resizing image
    image_.resize(camera_.width_, camera_.height_);

    /** \todo
     * (optional) Implement supersampling to avoid aliasing artifacts:
     * - Instead of casting just one, you cast n*n rays per pixel and store the average of the traced color values.
     * - Start with 2x2 rays per pixel and check the result with the cube scene.
     * - Try to generalize this to support arbitrary n*n supersampling
     * - To cast this many rays will slow down your image computation. To avoid this, implement adaptive supersampling:
     *      * Supersampling is just needed if neighboring pixels have a noticable difference in their color values
     *      * Start with one ray per pixel and store the result in a temporary `Image` variable
     *      * Now, iterate a second time over your image and compare the color values of the each pixel with those of its neighbors' pixels
     *      * If this difference is higher than a certain threshold, you apply 4x4 supersampling (SSAA16) for this pixel
     *      * Experiment with the rings scene to find a good threshold
     * Hints:
     * - Some image viewers will blur your images by default to avoid aliasing if you zoom in,
     * open your images with an image manipulation tool instead.
     * - It may help to visualize the x and y coordinates of your subrays on a sheet of paper
     * - Try to avoid sqrt computations like in norm when you want to compute a color difference, use normSq instead
     */

    // uncomment the following line to use multithreading
    #pragma omp parallel for
    for (int x = 0; x < camera_.width_; ++x)
    {
        for (int y = 0; y < camera_.height_; ++y)
        {
            // Ray ray = camera_.primary_ray(x, y);
            // vec3 color = trace(ray, 0);

            vec3 color(0.0);
            vec3 minColor(1.0), maxColor(0.0);

            double offsets[2][2][2];
            genOffsets<2>(offsets);

            vec3 colors[2][2];

            for(size_t i = 0; i < 2; i++) {
                for(size_t j = 0; j < 2; j++) {
                    Ray ray = camera_.primary_ray(x + offsets[i][j][0],
                                                    y + offsets[i][j][1]);

                    colors[i][j] = trace(ray, 0);
                    color += colors[i][j] / 4;
                    minColor = min(minColor, colors[i][j]);
                    maxColor = max(maxColor, colors[i][j]);
                }
            }

            bool repeatSS = norm(maxColor - minColor) > 0.02;

            if(repeatSS) {
                color = vec3(0.0);

                double offsets[4][4][2];
                genOffsets<4>(offsets);

                vec3 colors[4][4];

                for(size_t i = 0; i < 4; i++) {
                    for(size_t j = 0; j < 4; j++) {
                        Ray ray = camera_.primary_ray(x + offsets[i][j][0],
                                                        y + offsets[i][j][1]);

                        colors[i][j] = trace(ray, 0);
                        color += colors[i][j] / 16;
                    }
                }
            }

            // avoid over-saturation
            color = min(color, vec3(1, 1, 1));

            // store pixel color
            image_(x, y) = color;
        }
    }

}

//-----------------------------------------------------------------------------

void Raytracer::write_image(const std::string& filename)
{
    image_.write_png(filename);
}

//-----------------------------------------------------------------------------

vec3 Raytracer::trace(const Ray& ray, int depth)
{
    // stop if recursion depth (=number of reflection) is too large
    if (depth > max_depth_)
        return vec3(0, 0, 0);

    // Find first intersection with an object. If an intersection is found,
    // it is stored in object, point, normal, and t.
    Material material;
    vec3 point;
    vec3 normal;
    double t;
    if (!intersect_scene(ray, material, point, normal, t))
    {
        return background_;
    }

    // compute local Phong lighting (ambient+diffuse+specular)

    vec3 color = lighting(point, normal, -ray.direction_, material);

    if(material.mirror > 0.0) {
        vec3 reflectedDirection = mirror(-ray.direction_, normal);
        Ray reflectedRay(point, reflectedDirection);
        vec3 reflectedColor = trace(reflectedRay, depth + 1);

        color = (1 - material.mirror) * color + material.mirror * reflectedColor;
    }


    return color;
}

//-----------------------------------------------------------------------------

bool Raytracer::intersect_scene(const Ray& ray, Material& intersection_material,
                                vec3& intersection_point,
                                vec3& intersection_normal,
                                double& intersection_distance)
{
    double t, tmin(DBL_MAX);
    vec3 p, n, d;

    for (Object_ptr o : objects_) // for each object
    {
        if (o->intersect(ray, p, n, d, t)) // does ray intersect object?
        {
            if (t < tmin) // is intersection point the currently closest one?
            {
                tmin = t;
                intersection_material = o->material_;
                intersection_material.diffuse = d;
                intersection_point = p;
                intersection_normal = n;
                intersection_distance = t;
            }
        }
    }

    return (tmin < DBL_MAX);
}

//-----------------------------------------------------------------------------

bool Raytracer::shadow(const vec3& point, const vec3& light) {
    Ray shadowRay(point, light - point);
    Material temp;
    vec3 intersection;
    vec3 intersectionNormal;
    double t = -1.0;
    double distanceSq = normSq(point - light);

    intersect_scene(shadowRay, temp, intersection, intersectionNormal, t);

    return t > 0.0 && distanceSq > t * t;
}

vec3 Raytracer::lighting(const vec3& point, const vec3& normal,
                         const vec3& view, const Material& material)
{
    vec3 color(0.0, 0.0, 0.0);
 
    color += ambience_ * material.ambient;

    for(Light light: lights_) {
        if(material.shadowable && shadow(point, light.position))
            continue;

        vec3 l = normalize(light.position - point);
        double cosTheta = dot(normal, l);

        vec3 diffColor(0.0);
        vec3 specColor(0.0);

        if(cosTheta >= 0.0) {
            diffColor += material.diffuse * cosTheta;

            vec3 r = normalize(mirror(l, normal));
            double cosAlpha = dot(r, view);

            if(cosAlpha >= 0.0)
                specColor += material.specular
                            * std::pow(cosAlpha, material.shininess);
        }

        color += light.color * (diffColor + specColor);
    }

    return color;
}

//=============================================================================
