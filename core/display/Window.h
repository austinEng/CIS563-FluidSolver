//
// Created by austin on 2/25/16.
//

#ifndef FLUID_SIMULATOR_WINDOW_H
#define FLUID_SIMULATOR_WINDOW_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "InputHandler.h"

class Window {
public:
    Window(const char* title);
    Window(int w = 640, int h = 480);
    Window(int w, int h, const char* title);
    ~Window();
    void initloop();

private:
    GLFWwindow* _window;
    void setupInputCBs();
    void handleMouseInput(InputHandler::MouseState &mouseState);
};


#endif //FLUID_SIMULATOR_WINDOW_H
