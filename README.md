Ray Tracer
==========

This package contains exercise code for the lecture "Computer Graphics" by Prof. Dr. Mario Botsch, TU Dortmund.

Documentation
-------------

A pre-built HTML documentation can be found in `doc/index.html` and can be opened via any web browser.


Prerequisites
-------------

We use [CMake](www.cmake.org) for setting up build environments. You can download it [here](https://cmake.org/download/) or by using your favourite package manager.
You will also need a C++ compiler. Depending on your OS, the typical choices are:

  * Linux - GCC/G++, usually pre-installed (terminal command: `sudo apt install build-essential`)
  * MacOS - Clang, comes with the IDE "XCode" (terminal command: `xcode-select --install`)
  * Windows - MSVC, can be installed alongside the IDE "Visual Studio Community" (see below) or separately [here](https://visualstudio.microsoft.com/downloads/#other) under "All Downloads" --> "Tools for Visual Studio" --> "Build Tools for Visual Studio 2019"

We highly recommand to use some kind of IDE. Prominent examples are:
 
 - [XCode](https://developer.apple.com/xcode/) (MacOS) 
 - [Visual Studio Community](https://visualstudio.microsoft.com/de/vs/community/) (Windows)
 - [Visual Studio Code](https://code.visualstudio.com/) (Linux, Windows, MacOS)
 - [Jetbrains CLion](https://www.jetbrains.com/de-de/clion/) (Linux, Windows, MacOS)
 
Below, we provide examples for setting up, compiling and running the project via command line, XCode, VSCommunity and VSCode.

Building on Linux via Command Line (no IDE)
-------------------------------------------

Execute the following commands in the exercise's top-level directory:

    mkdir build
    cd build
    cmake ..
    make

The last command (`make`) compiles the application. Re-run it whenever you have added/changed code in order to re-compile.

For running the code via command line use

    ./raytracer


Building on MacOS (XCode)
--------------------------

Execute the following commands in the exercise's top-level directory:

    mkdir xcode
    cd xcode
    cmake -G Xcode ..

Open the generated .xcodeproj file. Inside XCode, select the raytracer-executable in the top bar next to the stop button. You can run the code by pressing the play button.
You can specify command line arguments by again opening the menu next to the stop button from the top bar and selecting "Edit Scheme". In the popup window, select the "Run" option left and open the "Arguments" tab.


Building on Windows (Visual Studio Community) 
---------------------------------------------

* In order to get the MSVC C++ compiler, make sure that you check "Desktop development with C++" during installation of [VSCommunity](https://visualstudio.microsoft.com/de/vs/community/)
* Create an empty build folder inside the project's top-level directory
* Start cmake-gui.exe (located in cmake's bin folder)
* Specify the top-level directory as source directory (button Browse source...)
* Specify the previously created build folder as build directory (button Browse build...)
* Select "Configure" using your Visual Studio Version as option.
* When configuration is finished, select "Generate".
* Start Visual Studio Community
* Open the project via File -> open -> project -> .sln in build folder
* In the project explorer window on the right, right-click the project (raytracer) and set it as startup-project
* Switch to release mode
* You can specify command line arguments via project -> properties -> debugging -> command arguments
* Hit CTRL + F5 to build and run (or CTRL + SHIFT + B to build)


Building via VSCode
-------------------

There are a lot of useful extensions for VSCode to shape it the way you like.
The required extensions for C++ development and Cmake support are "C/C++" and "Cmake Tools". Extensions can be found in the extensions tab on the left.
Once this is done, you can set up the project:

 * Start VSCode
 * Open the project via File --> Open Folder and select the exercise's top-level directory containing the `CMakeLists.txt` file and accept by clicking `OK`
 * In the bottom bar, click on `CMake`, choose your compiler and select `Release` mode
 * In the bottom bar, click on  `[all]` and select the right target (raytracer)
 * Still in the bottom bar, there are buttons for building (Build) and running (play symbol)
 * You can specify command line arguments by using the terminal inside of VSCode

Have a look at the VSCode [documentation](https://code.visualstudio.com/docs/cpp/introvideos-cpp) for further details. 


Command Line Arguments
----------------------

Running the ray tracer without any arguments will render all scenes. 
If you want to ray trace a specific scene, you can pass two command line arguments: an input scene (`*.sce`) and an output image (`*.png`).
For example, when using the command line, rendering the sphere scene can be invoked via

    ./raytracer ../scenes/spheres/spheres.sce output.png
    
Each IDE provides a different way, for specifying command line arguments. We give some examples above but refer to the IDEs' documentations for further details.



Code Overview
-------------

The best way to become familiar with the code is to look at the Raytracer class implemented in `Raytracer.h` and `Raytracer.cpp`. The function `load_scene()` loads all objects, lights and camera settings from a `.sce` file.  The function `compute_image()` is the starting point for the ray tracer that creates a ray for each pixel of the image and calls `trace()` to determine each pixel's color.


Recommended order for the exercises
-----------------------------------

We recommend that you work on the tasks in the following order. All tasks can be found in the ToDo list as well.

Note that, to give you something to compare to, we provide output images created with a completed version of the exercise. They can be found in the respective scene directories, e.g. `scenes/spheres/`. Refer to the lecture slides to get an idea of what the intermediate steps look like.

* **Phong Shading and Reflections:**
Begin by implementing the Phong shading model. Your starting point will be the `lighting()` function in `Raytracer.cpp`. Once you have finished the Phong shading and get a realistic illumination, continue with reflections, which have to be added to `trace()` in `Raytracer.cpp`. During coding, ensure you are on the right track by ray tracing the spheres (`scenes/spheres/spheres.sce`).

* **Intersection computation with a plane:**
To complete proper rendering of the spheres, implement ray-plane intersections (`Plane::intersect()` in `Plane.cpp`).

* **Triangle Meshes:**
Ray tracing spheres is great (try the scenes with molecules!), but because we want to be able to render more complex scenes, we cannot stop there. The next step is the ray-triangle intersection (`Mesh::intersect_triangle()` in `Mesh.cpp`). Start testing with the cube (`scenes/cube/cube.sce`), which uses flat shading. Then try the office (`scenes/office/office.sce`) or the toon faces (`scenes/toon_faces/toon_faces.sce`). For those scenes, computation of vertex normals needs to be implemented (`Mesh::compute_normals` in `Mesh.cpp`). Implement weighted vertex normals to get a nice result in the rings scene (`scenes/rings/rings.sce`).

* **Acceleration:**
Are you annoyed because ray-tracing is damn slow? Let's make things faster! First implement the bounding box test for triangle meshes (`Mesh::intersect_bounding_box()` in `Mesh.cpp`). With this feature you should see a significant performance boost rendering the toon faces (`scenes\toon_faces\toon_faces.sce`). Also try using OpenMP to enable multi-core parallelization. For this, simply put
	(`#pragma omp parallel for`)
above the outer loop in during ray generation (`compute_image()` in `Raytracer.cpp`).

* **Textures (optional)**
One color per object is boring. Complete the texture support in `Mesh::intersect_triangle()` and render our cool Pokemon scene (`scenes\pokemon\pokemon.sce`).

* **Supersampling (optional):**
Do you still want to explore more features? Why not try super-sampling in order to reduce the aliasing problems? Simply adjust the function `compute_image()`, such that you shoot several rays for each pixel and average their resulting colors. To avoid loosing too much performance, you should implement adaptive supersampling (details in `compute_image()`). There are also 16-ray-supersampled versions of the images in each scene folder for comparison (suffix _SS16). 

* **Render all scenes:** 
If you finished all the tasks above, you can render all scenes (no command line argument).
Compare your results to those in the scenes folder. You should get the same pictures, if you did no mistakes.


License
-------

The code examples are copyright Computer Graphics Group, TU Dortmund.


Have fun!

