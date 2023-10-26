file(GLOB_RECURSE SRCS ./*.cpp)
file(GLOB_RECURSE HDRS ./*.h)

add_executable(raytracer ${SRCS} ${HDRS} )
target_link_libraries(raytracer stb_image)
