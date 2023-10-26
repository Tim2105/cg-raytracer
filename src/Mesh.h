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
#include "utils/Image.h"
#include <vector>
#include <float.h>
#include <string>

//== CLASS DEFINITION =========================================================

/// \class Mesh Mesh.h
/// This class represents a simple triangle mesh, stored as an indexed face set,
/// i.e., as an array of vertices and an array of triangles.
class Mesh : public Object
{
public:
    /// This type is used to choose between flat shading and Phong shading
    enum Draw_mode
    {
        FLAT,
        PHONG
    };

    /// Construct a mesh by specifying its draw_mode and its filename
    Mesh(Draw_mode draw_mode = FLAT, std::string filename = 0);

    /// Intersect mesh with ray (calls ray-triangle intersection)
    /// If \c ray intersects a face of the mesh, it provides the following results:
    /// \param[in] ray the ray to intersect the mesh with
    /// \param[out] intersection_point the point of intersection
    /// \param[out] intersection_normal the surface normal at intersection point
    /// \param[out] intersection_distance ray parameter at the intersection point
    virtual bool intersect(const Ray& ray, vec3& intersection_point,
                           vec3& intersection_normal,
                           vec3& intersection_diffuse,
                           double& intersection_distance) const;

private:
    /// a vertex consists of a position and a normal
    struct Vertex
    {
        /// vertex position
        vec3 position;
        /// vertex normal
        vec3 normal;
    };

    /// a triangle is specified by three indices and a normal
    struct Triangle
    {
        /// index of first vertex (for array Mesh::vertices_)
        int i0;
        /// index of second vertex (for array Mesh::vertices_)
        int i1;
        /// index of third vertex (for array Mesh::vertices_)
        int i2;
        /// triangle normal
        vec3 normal;

        /// texture-coordinate (uv) index of first, second and third vertex (for array Mesh::vertices_)
        int iuv0;
        int iuv1;
        int iuv2;
    };

private:
    /// Read mesh form an OBJ file
    bool read_obj(const char* filename);

    /// Read mtl file of an obj object
    bool read_mtl(std::string path, std::vector<Image>& textures);

    /// Compute normal vectors for triangles and vertices
    void compute_normals();

    /// Compute the axis-aligned bounding box, store minimum and maximum point in bb_min_ and bb_max_
    void compute_bounding_box();

    /// Does \c _ray intersect the bounding box of the mesh?
    bool intersect_bounding_box(const Ray& ray) const;

    /// Intersect a triangle with a ray. Return whether there is an intersection.
    /// If there is an intersection, store intersection data.
    /// This function overrides Object::intersect().
    /// \param[in] triangle the triangle to be intersected
    /// \param[in] ray the ray to intersect the triangle with
    /// \param[out] intersection_point the point of intersection
    /// \param[out] intersection_normal the surface normal at intersection point
    /// \param[out] intersection_diffuse the diffuse color at intersection point
    /// \param[out] intersection_distance ray parameter at the intersection point
    bool intersect_triangle(const Triangle& triangle, const Ray& ray,
                            vec3& intersection_point, vec3& intersection_normal,
                            vec3& intersection_diffuse,
                            double& intersection_distance) const;

public:
    /// Does this mesh use flat or Phong shading?
    Draw_mode draw_mode_;

    /// Array of vertices
    std::vector<Vertex> vertices_;

    /// Array of triangles
    std::vector<Triangle> triangles_;

    /// Minimum point of the bounding box
    vec3 bb_min_;
    /// Maximum point of the bounding box
    vec3 bb_max_;

    /// mesh's texture and uv-coordinates
    Image texture_;

    /// u-coordinates to access texture from 0 (left) to 1 (right)
    std::vector<double> u_coordinates_;

    /// v-coordinates to access texture from 0 (bottom) to 1 (top)
    std::vector<double> v_coordinates_;

    /// indicates if mesh has a texture
    bool hasTexture_;
};

//-----------------------------------------------------------------------------

/// read sphere from stream
inline std::istream& operator>>(std::istream& is, Mesh& m)
{
    std::string filename, mode;
    is >> filename >> mode;

    m = Mesh((mode == "FLAT" ? Mesh::FLAT : Mesh::PHONG), filename.c_str());
    is >> m.material_;

    return is;
}

//=============================================================================
