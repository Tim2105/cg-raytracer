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

        vec3 p0ToP1 = p1 - p0;
        vec3 p0ToP2 = p2 - p0;
        vec3 p1ToP2 = p2 - p1;

        vertices_[t.i0].normal += angle(p0ToP1, p0ToP2) * t.normal;
        vertices_[t.i1].normal += angle(p0ToP1, p1ToP2) * t.normal;
        vertices_[t.i2].normal += angle(p0ToP2, p1ToP2) * t.normal;
    }

    // Normiere alle Vertex-Normalen
    for (Vertex& v : vertices_)
        v.normal = normalize(v.normal);

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
    double t1, t2, tMin = 0.0, tMax = DBL_MAX;

    // Durchlaufe alle Koordinaten
    for(size_t i = 0; i < 3; i++) {
        // Sonderfall: d_i = 0
        if(std::abs(ray.direction_[i]) <= DBL_MIN) {
            if(ray.origin_[i] < bb_min_[i] || ray.origin_[i] > bb_max_[i])
                return false;
        } else {
            // Schnittdistanz 1
            t1 = (bb_min_[i] - ray.origin_[i]) / ray.direction_[i];
            // Schnittdistanz 2
            t2 = (bb_max_[i] - ray.origin_[i]) / ray.direction_[i];

            // Sorge für t1 <= t2
            if(t1 > t2)
                std::swap(t1, t2);

            // Überschreibe t1 in tMin, wenn t1 > 0
            if(t1 > tMin)
                tMin = t1;
            
            tMax = t2;

            // Vergleiche tMin und tMax (müssen auch >= 0 sein)
            if(tMin > tMax)
                return false;
            if(tMax < 0.0)
                return false;
        }
    }

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
    intersection_distance = t;

    // Berechne die Normale des Schnittpunkts
    if(draw_mode_ == PHONG) {
        // Smooth Shading
        double alpha = 1.0 - beta - gamma;

        vec3 interpolNormal(0.0);
        interpolNormal += alpha * vertices_[triangle.i0].normal;
        interpolNormal += beta * vertices_[triangle.i1].normal;
        interpolNormal += gamma * vertices_[triangle.i2].normal;

        intersection_normal = normalize(interpolNormal);
    } else // Flat Shading
        intersection_normal = triangle.normal;

    if(hasTexture_) {
        double alpha = 1.0 - beta - gamma;
        double uRel = alpha * u_coordinates_[triangle.iuv0] + 
                        beta * u_coordinates_[triangle.iuv1] +
                        gamma * u_coordinates_[triangle.iuv2];

        double vRel = alpha * v_coordinates_[triangle.iuv0] + 
                        beta * v_coordinates_[triangle.iuv1] +
                        gamma * v_coordinates_[triangle.iuv2];

        int uPixel = uRel * texture_.width();
        int vPixel = vRel * texture_.height();

        intersection_diffuse = texture_(uPixel, vPixel);
    } else
        intersection_diffuse = material_.diffuse;

    return true;
}

//=============================================================================
