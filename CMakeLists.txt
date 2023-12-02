cmake_minimum_required(VERSION 3.5)
project(nifs3edit C)

set(OpenGL_GL_PREFERENCE LEGACY)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

add_executable(nifs3edit nifs3edit.c)
target_link_libraries(nifs3edit ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} m)
target_include_directories(nifs3edit PUBLIC ${OPENGL_INCLUDE_DIRS}
                                            ${GLUT_INCLUDE_DIRS} .)
