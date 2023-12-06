cmake_minimum_required(VERSION 3.5)
project(nifs3edit C)

set(OpenGL_GL_PREFERENCE LEGACY)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

add_executable(nifs3edit nifs3edit.c)
target_link_libraries(nifs3edit OpenGL::GL GLUT::GLUT m)
target_include_directories(nifs3edit PUBLIC .)
