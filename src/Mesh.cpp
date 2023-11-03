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

#include "Mesh.h"
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <map>

//== IMPLEMENTATION ===========================================================

Mesh::Mesh(Draw_mode _draw_mode, std::string _filename)
{
    // set draw mode
    draw_mode_ = _draw_mode;

    hasTexture_ = false;

    // load mesh from file
    read_obj(_filename.c_str());
}

bool Mesh::read_obj(const char* _filename)
{
    // open obj file
    std::ifstream ifs(_filename);
    if (!ifs)
    {
        std::cerr << "Can't open " << _filename << "\n";
        return false;
    }

    bool hasNormals = false;
    bool hasUV = false;

    std::string filename(_filename);
    std::string line;
    int counter = -1;
    std::map<int, bool> uvDone;
    std::vector<Image> textures;
    // parse line by line
    while (std::getline(ifs, line))
    {
        //devide line into header (first word) and lineData (rest)
        size_t firstSpace = line.find_first_of(" ");
        std::string header = line.substr(0, firstSpace);
        std::istringstream lineData(line.substr(firstSpace + 1));

        //vertices
        if (header == "v")
        {
            Vertex v;
            lineData >> v.position[0] >> v.position[1] >> v.position[2];

            vertices_.push_back(v);
            continue;
        }

        //uv-coordinates
        if (header == "vt")
        {
            hasUV = true;

            double u, v;

            lineData >> u >> v;

            if (u > 1.0 || u < 0.0)
                u -= floor(u);
            if (v > 1.0 || v < -0.0)
                v -= floor(v);

            u_coordinates_.push_back(u);
            v_coordinates_.push_back(v);
            continue;
        }

        if (header == "vn")
        {
            hasNormals = true;
            continue;
        }

        // material file
        if (header == "mtllib")
        {
            std::stringstream mtl;
            mtl << filename.substr(0, filename.find_last_of("/") + 1)
                << lineData.str();

            if (!read_mtl(mtl.str(), textures))
            {
                std::cerr << "Cannot read mtl file " << mtl.str() << std::endl;
            }

            if (textures.size() > 0)
                hasTexture_ = true;

            continue;
        }

        // start of new material
        if (header == "usemtl")
        {
            counter++;
            continue;
        }

        // faces
        if (header == "f")
        {
            Triangle t;

            int uv[3];

            enum
            {
                NORMALS,
                UV,
                BOTH,
                NONE
            } nuv_status;
            if (hasUV)
                nuv_status = hasNormals ? BOTH : UV;
            else
                nuv_status = hasNormals ? NORMALS : NONE;

            // dummy varaibles for / and normal indices
            int d1;
            char d2;

            // read in face indices and uv indices, skip normal indices
            switch (nuv_status)
            {
                case BOTH:
                    // format: index0/texture0/normal0 index1/texture1/normal1 index2/texture2/normal2
                    lineData >> t.i0 >> d2 >> uv[0] >> d2 >> d1;
                    lineData >> t.i1 >> d2 >> uv[1] >> d2 >> d1;
                    lineData >> t.i2 >> d2 >> uv[2] >> d2 >> d1;

                    uv[0]--;
                    uv[1]--;
                    uv[2]--;
                    t.iuv0 = uv[0];
                    t.iuv1 = uv[1];
                    t.iuv2 = uv[2];
                    break;
                case UV:
                    // format: index0/texture0 index1/texture1 index2/texture2
                    lineData >> t.i0 >> d2 >> uv[0];
                    lineData >> t.i1 >> d2 >> uv[1];
                    lineData >> t.i2 >> d2 >> uv[2];

                    uv[0]--;
                    uv[1]--;
                    uv[2]--;
                    t.iuv0 = uv[0];
                    t.iuv1 = uv[1];
                    t.iuv2 = uv[2];
                    [[fallthrough]];
                case NORMALS:
                    // format: index0//normal0 index1//normal1 index2//normal2
                    lineData >> t.i0 >> d2 >> d2 >> d1;
                    lineData >> t.i1 >> d2 >> d2 >> d1;
                    lineData >> t.i2 >> d2 >> d2 >> d1;
                    [[fallthrough]];
                case NONE:
                    // format: index0 index1 index2
                    lineData >> t.i0 >> t.i1 >> t.i2;
            }

            //decrement because obj indices start by 1
            t.i0--;
            t.i1--;
            t.i2--;

            //convert uv coordinates s.th. we can use just one big combined tex instead of multiple ones
            for (int i = 0; i < 3 && hasUV; i++)
            {
                if (!uvDone[uv[i]])
                {
                    int combinedW = 0;
                    for (int i = 0; i < counter; i++)
                    {
                        combinedW += textures[i].width();
                    }
                    u_coordinates_[uv[i]] =
                        (u_coordinates_[uv[i]] * textures[counter].width() +
                         combinedW) /
                        static_cast<double>(texture_.width());
                    v_coordinates_[uv[i]] =
                        (v_coordinates_[uv[i]] * textures[counter].height()) /
                        static_cast<double>(texture_.height());
                    uvDone[uv[i]] = true;
                }
            }

            // add triangle to our triangle vector
            triangles_.push_back(t);
        }
    }

    std::cout << "\n  read " << _filename << ": " << vertices_.size()
              << " vertices, " << triangles_.size() << " triangles"
              << std::flush;

    compute_bounding_box();
    compute_normals();

    return true;
}

//-----------------------------------------------------------------------------

bool Mesh::read_mtl(std::string path, std::vector<Image>& textures)
{
    // open mtl file
    std::ifstream ifs(path.c_str());
    if (!ifs)
    {
        std::cerr << "Can't open " << path << "\n";
        return false;
    }

    std::string line;
    int numTexturesPerMaterial = 0;
    Image tmpimage;

    // parse line by line
    while (std::getline(ifs, line))
    {
        //devide line into header (first word) and lineData (rest)
        size_t firstSpace = line.find_first_of(" ");
        std::string header = line.substr(0, firstSpace);
        std::istringstream lineData(line.substr(firstSpace + 1));

        if (header == "newmtl")
        {
            numTexturesPerMaterial = 0;
            continue;
        }

        if (header.substr(0, 3) == "map" && numTexturesPerMaterial == 0)
        {
            std::stringstream tmp;
            tmp << path.substr(0, path.find_last_of("/") + 1) << lineData.str();

            tmpimage.read(tmp.str().c_str());
            textures.push_back(tmpimage);
            numTexturesPerMaterial++;
        }
    }

    unsigned int maxH = 0;
    unsigned int sumW = 0;
    for (size_t i = 0; i < textures.size(); i++)
    {
        sumW += textures[i].width();
        maxH = std::max(maxH, textures[i].height());
    }
    texture_.resize(sumW, maxH);

    for (unsigned int x = 0; x < sumW; x++)
    {
        for (unsigned int y = 0; y < maxH; y++)
        {
            unsigned int texnr = 0;
            unsigned int combinedW = 0;

            for (size_t i = 0; i < textures.size(); i++)
            {
                if (x >= combinedW + textures[i].width())
                {
                    combinedW += textures[i].width();
                    texnr++;
                }
                else
                    break;
            }

            if (y < textures[texnr].height())
            {
                texture_(x, y) = textures[texnr](x - combinedW, y);
            }
            else
            {
                texture_(x, y) = vec3(0, 0, 0);
            }
        }
    }

    return true;
}

//-----------------------------------------------------------------------------

void Mesh::compute_normals()
{
    /** \todo
     * In some scenes (e.g the office scene) some objects have to be flat
     * shaded (e.g. the desk) while other objects should be Phong shaded to appear
     * realistic (e.g. chairs). We provide code to compute triangle normals. You have to
     * implement the following:
     * - Compute vertex normals by averaging the normals of their incident triangles.
     * - Store the vertex normals in the Vertex::normal member variable.
     * - Weigh the normals by their triangles' angles.
     */


    // initialize vertex normals to zero
    for (Vertex& v : vertices_)
    {
        v.normal = vec3(0, 0, 0);
    }

    // compute triangle normals
    for (Triangle& t : triangles_)
    {
        const vec3& p0 = vertices_[t.i0].position;
        const vec3& p1 = vertices_[t.i1].position;
        const vec3& p2 = vertices_[t.i2].position;
        t.normal = normalize(cross(p1 - p0, p2 - p0));
    }

}

//-----------------------------------------------------------------------------

void Mesh::compute_bounding_box()
{
    bb_min_ = vec3(DBL_MAX, DBL_MAX, DBL_MAX);
    bb_max_ = vec3(-DBL_MAX, -DBL_MAX, -DBL_MAX);

    for (Vertex v : vertices_)
    {
        bb_min_ = min(bb_min_, v.position);
        bb_max_ = max(bb_max_, v.position);
    }
}

//-----------------------------------------------------------------------------

bool Mesh::intersect_bounding_box(const Ray& ray) const
{
    UNUSED(ray);

    /** \todo
    * Intersect the ray `_ray` with the axis-aligned bounding box of the mesh.
    * Note that the minimum and maximum point of the bounding box are stored
    * in the member variables `bb_min_` and `bb_max_`. Return whether the ray
    * intersects the bounding box.
    *
    * Hints:
    * - The box intersection is basically a combination of 6 plane/ray intersections.
    * - The resulting intersection point has to be checked against the lower and upper bounds `bb_min_` and `bb_max_`.
    * - The intersection tests can be done easier than the existing general ray/plane intersection because the planes are axis aligned.
    * - One positive intersection is sufficient to return true.
    * - It helps to visualize a 2D ray/rectangle intersection on a sheet of paper.
    * - To debug your bounding box code, comment out your triangle intersection test (return true at the beginning) and test
    *   the cube or toon_faces scene. You should see the black bounding boxes, if everything is correct.
    * - There are other (faster) ray/box intersection approaches, feel free to implement one of those instead.
    *
    * Note:
    * This function is used in `Mesh::intersect()` to avoid the intersection test
    * with all triangles of every mesh in the scene. The bounding boxes are computed
    * in `Mesh::compute_bounding_box()`.
    */


    return true;
}

//-----------------------------------------------------------------------------

bool Mesh::intersect(const Ray& ray, vec3& intersection_point,
                     vec3& intersection_normal, vec3& intersection_diffuse,
                     double& intersection_distance) const
{
    // check bounding box intersection
    if (!intersect_bounding_box(ray))
    {
        return false;
    }

    vec3 p, n, d;
    double t = 0.0;

    intersection_distance = DBL_MAX;

    // for each triangle
    for (const Triangle& triangle : triangles_)
    {
        // does ray intersect triangle?
        if (intersect_triangle(triangle, ray, p, n, d, t))
        {
            // is intersection closer than previous intersections?
            if (t < intersection_distance)
            {
                // store data of this intersection
                intersection_distance = t;
                intersection_point = p;
                intersection_normal = n;
                intersection_diffuse = d;
            }
        }
    }

    return (intersection_distance != DBL_MAX);
}

//-----------------------------------------------------------------------------
class Mat3x3 {
    private:
        double data[9] = {0};

    public:
        constexpr Mat3x3(const vec3& c1, const vec3& c2, const vec3& c3) {
            data[0] = c1[0];
            data[1] = c2[0];
            data[2] = c3[0];
            data[3] = c1[1];
            data[4] = c2[1];
            data[5] = c3[1];
            data[6] = c1[2];
            data[7] = c2[2];
            data[8] = c3[2];
        }

        /**
         * @brief Berechne die Determinante mit der Regel von Sarrus
         */
        constexpr double determinant() {
            return data[0] * data[4] * data[8] +
                   data[1] * data[5] * data[6] +
                   data[2] * data[3] * data[7] -
                   data[2] * data[4] * data[6] -
                   data[1] * data[3] * data[8] -
                   data[0] * data[5] * data[7];
        }
};

bool Mesh::intersect_triangle(const Triangle& triangle, const Ray& ray,
                              vec3& intersection_point,
                              vec3& intersection_normal,
                              vec3& intersection_diffuse,
                              double& intersection_distance) const
{
    intersection_diffuse = material_.diffuse;
    
    // Berechne die Lösungen des Systems
    // origin - pA = beta(pB - pA) + gamma(pC - pA) - t * direction
    vec3 pA = vertices_[triangle.i0].position;
    vec3 pB = vertices_[triangle.i1].position;
    vec3 pC = vertices_[triangle.i2].position;

    // Mithilfe der Cramer'schen Regel
    Mat3x3 mA(pB - pA, pC - pA, -ray.direction_);
    Mat3x3 mA1(ray.origin_ - pA, pC - pA, -ray.direction_);
    Mat3x3 mA2(pB - pA, ray.origin_ - pA, -ray.direction_);
    Mat3x3 mA3(pB - pA, pC - pA, ray.origin_ - pA);

    double detA = mA.determinant();
    
    // Wenn die Determinante von mA = 0 ist, hat das LGS keine Lösung
    if(detA == 0.0)
        return false;

    double beta = mA1.determinant() / detA;
    double gamma = mA2.determinant() / detA;
    double t = mA3.determinant() / detA;

    // Wenn beta oder gamma negativ sind, oder beta + gamma > 1,
    // liegt der Schnittpunkt nicht im Dreieck
    if(beta < 0.0 || gamma < 0.0 || beta + gamma > 1.0)
        return false;

    // Wenn t < Eplison, schneidet der Strahl nicht in positive Richtung
    if(t <= 1e-5)
        return false;

    intersection_point = ray(t);
    intersection_normal = triangle.normal;
    intersection_distance = t;

    return true;

    /** \todo
    * Intersect ray with triangle:
    * - store intersection point in `intersection_point`
    * - store ray parameter in `intersection_distance`
    * - store normal at intersection point in `intersection_normal`.
    * - Depending on the member variable `draw_mode_`, use either the triangle
    *  normal (`Triangle::normal`) or interpolate the vertex normals (`Vertex::normal`).
    * - return `true` if there is an intersection with t > 1e-5 (avoid "shadow acne")
    *
    * Hint: Rearrange `ray.origin + t*ray.dir = a*p0 + b*p1 + (1-a-b)*p2` to obtain a solvable
    * system for a, b and t.
    * Refer to [Cramer's Rule](https://en.wikipedia.org/wiki/Cramer%27s_rule) to easily solve it.
     */

    // const vec3& p0 = vertices_[triangle.i0].position;
    // const vec3& p1 = vertices_[triangle.i1].position;
    // const vec3& p2 = vertices_[triangle.i2].position;

        /** \todo
    * (optional) Support textured triangles:
    * - `hasTexture_` indicates if the mesh is textured.
    * - Acces the three u and v texture coordinates stored in `u_coordinates` resp. `v_coordinates` via the triangles iuv indices.
    * - Interpolate the uv-coordinates using your barycentric coordinates to get the intersection point's uv.
    * - Convert the relative uv coordinates (from 0 to 1) to absolute pixel coordinates (from 0 to width/height - 1 of `texture_`)
    * - Store the resulting texture color in `intersection_diffuse`
    * - You will notice that there will be shadows on the sky mesh in the pokemon scene.
    * Use `material.shadowable` in the `lighting(...)` function to prevent it from being shadowed.
    * (`material.shadowable` is already set to false for the sky mesh and true for all other meshes, so you don't have to set it by yourself)
     */
}

//=============================================================================
