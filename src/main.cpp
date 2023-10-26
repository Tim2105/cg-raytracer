//=============================================================================
//
//   Exercise code for the lecture
//   "Computer Graphics"
//   by Prof. Dr. Mario Botsch, TU Dortmund
//
//   Copyright (C) Computer Graphics Group, TU Dortmund.
//
//=============================================================================

//== defines =================================================================

#ifndef SCENE_PATH
#if defined(_WIN32) || defined(__XCODE__)
#define SCENE_PATH ("../../scenes/")
#else
#define SCENE_PATH ("../scenes/")
#endif
#endif

//== includes =================================================================

#include "Raytracer.h"
#include "utils/StopWatch.h"

int main(int argc, char** argv)
{
    if (argc == 2 || argc > 3)
    {
        std::cerr << "Usage: " << argv[0] << "  input.sce output.png\n\n"
                  << std::endl;
        exit(1);
    }

    if (argc == 3)
    {
        // read scene from file
        std::cout << "Read scene...";
        Raytracer raytracer(argv[1]);

        // raytrace image
        StopWatch timer;
        timer.start();
        std::cout << "Ray tracing..." << std::flush;

        raytracer.compute_image();
        timer.stop();
        std::cout << "done (" << timer << ")\n";

        // write the resulting image
        std::cout << "Write image...";
        raytracer.write_image(argv[2]);
        std::cout << "done\n";
    }

    if (argc == 1)
    {
        std::cout << "rendering all scenes..." << std::endl;

        std::string path = SCENE_PATH;
        for (int i = 0; i < 10; i++)
        {
            std::string si, so;
            switch (i)
            {
                case 0:
                    si = path + std::string("spheres/spheres.sce");
                    so = "o_01_spheres.png";
                    break;
                case 1:
                    si = path + std::string("shadow/shadow.sce");
                    so = "o_02_shadow.png";
                    break;
                case 2:
                    si = path + std::string("mirror/mirror.sce");
                    so = "o_03_mirror.png";
                    break;
                case 3:
                    si = path + std::string("molecule/molecule.sce");
                    so = "o_04_molecule.png";
                    break;
                case 4:
                    si = path + std::string("cube/cube.sce");
                    so = "o_05_cube.png";
                    break;
                case 5:
                    si = path + std::string("mask/mask.sce");
                    so = "o_06_mask.png";
                    break;
                case 6:
                    si = path + std::string("toon_faces/toon_faces.sce");
                    so = "o_07_toon_faces.png";
                    break;
                case 7:
                    si = path + std::string("office/office.sce");
                    so = "o_08_office.png";
                    break;
                case 8:
                    si = path + std::string("rings/rings.sce");
                    so = "o_09_rings.png";
                    break;
                case 9:
                    si = path + std::string("pokemon/pokemon.sce");
                    so = "o_10_pokemon.png";
                    break;
            }

            // read scene from file
            std::cout << "Read scene..." << std::flush;
            Raytracer raytracer(si);

            // raytrace image
            StopWatch timer;
            timer.start();
            std::cout << "Ray tracing..." << std::flush;
            raytracer.compute_image();
            timer.stop();
            std::cout << "done (" << timer << ")\n";

            // write the resulting image
            std::cout << "Write image..." << std::flush;
            raytracer.write_image(so.c_str());
            std::cout << "done\n";
        }
    }
}
