# credit - base CMake config : Yining Karl Li , edited CMake config: Akshay Shah & Debanshu Singh

#name your project
project(FluidSolver)
cmake_minimum_required(VERSION 2.8)

# set creates a variable
set (NUPARU "nuparu")
# add include and src directories to path
include_directories(
        ${NUPARU}/include
        ${NUPARU}/src
)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Add path for pre-compiled libraries here (we will later link them with our compiled source)
# Add Nuparu library to path for OSX, linux and windows
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${NUPARU}/lib/osx)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${NUPARU}/lib/linux /usr/lib64)
elseif(WIN32)
    set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${NUPARU}/lib/win)
endif()

#set include variables
set(GLFW_INCLUDE_DIR ${NUPARU}/include )
set(GLEW_INCLUDE_DIR ${NUPARU}/include )
set(GLFW_LIBRARY_DIR ${CMAKE_LIBRARY_PATH})
set(GLEW_LIBRARY_DIR ${CMAKE_LIBRARY_PATH})

# Use find_package & find_library to link with
find_package(OpenGL REQUIRED)
find_package(GLEW)
find_library(GLFW_LIBRARY "glfw3" HINTS ${GLFW_LIBRARY_DIR})
find_library(JSONCPP "jsoncpp")

add_definitions(
        -DTW_STATIC
        -DTW_NO_LIB_PRAGMA
        -DTW_NO_DIRECT3D
        -DGLEW_STATIC
        -D_CRT_SECURE_NO_WARNINGS
)

set(CORE_LIBS ${GLFW_LIBRARY} ${GLUT_LIBRARY} ${GLEW_LIBRARY} ${JSONCPP} ${OPENGL_LIBRARY} )

# OSX-specific hacks/fixes
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    #Link IOKit because this is where we get GL stuff for OSX
    set(COCOA "-framework Cocoa")
    set(COREVIDEO "-framework CoreVideo")
    set(IOKIT "-framework IOKit")
    set(CORE_LIBS ${CORE_LIBS} ${COCOA} ${IOKIT} ${COREVIDEO})
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Linux specific hacks/fixes
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CORE_LIBS ${CORE_LIBS} X11 Xxf86vm Xrandr pthread Xi)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lX11 -lXxf86vm -lXrandr -lpthread -lXi")
endif()

# set compiler flags for c++11
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" OR ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -m64 -msse2 -w")
elseif(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

if(MSVC)
    set(COMPILER_FLAGS
            CMAKE_CXX_FLAGS
            CMAKE_CXX_FLAGS_DEBUG
            CMAKE_CXX_FLAGS_RELEASE
            CMAKE_C_FLAGS
            CMAKE_C_FLAGS_DEBUG
            CMAKE_C_FLAGS_RELEASE
        )
endif()

# Add source files you want to compile (.cpp)
set(CORE_SRC
        main.cpp
        core/display/InputHandler.cpp
        core/display/InputHandler.h
        core/display/Window.cpp
        core/display/Window.h
        core/fileIO/SceneLoader.cpp
        core/fileIO/SceneLoader.h
        core/geometry/Bound.cpp
        core/geometry/Bound.h
        core/geometry/Box.cpp
        core/geometry/Box.h
        core/geometry/Geo.cpp
        core/geometry/Geo.h
        core/geometry/GeoObject.h
        core/solver/FluidSolver.cpp
        core/solver/FluidSolver.h
        core/util/math.h
)

add_executable(FluidSolver ${CORE_SRC})
target_link_libraries(FluidSolver ${CORE_LIBS})
