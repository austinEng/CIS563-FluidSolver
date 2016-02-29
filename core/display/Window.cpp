//
// Created by austin on 2/25/16.
//

#include "Window.h"
#include <algorithm>
#include <iostream>

static void error_callback(int error, const char* description) {
    fputs(description, stderr);
}

InputHandler &inputHandler = InputHandler::getInputHandler();

Window::Window(const char *title) : Window(640, 480, title) { }
Window::Window(int w, int h) : Window(w, h, "GL Window"){ }
Window::Window(int w, int h, const char* title) : _window(nullptr) {
    glfwSetErrorCallback(error_callback);

    if (!glfwInit()) exit(EXIT_FAILURE);

    glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); //We don't want the old OpenGL

    _window = glfwCreateWindow(w, h, title, NULL, NULL);

    if (!_window) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(_window);
    glfwSwapInterval(1);
    setupInputCBs();

    glewExperimental= GL_TRUE; // Needed in core profile
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        exit(EXIT_FAILURE);
    }
}

Window::~Window() {

}

void Window::setupInputCBs() {

    glfwSetKeyCallback(_window, [&](GLFWwindow *window, int key, int scancode, int action, int mods) {
        std::cout << key << std::endl;
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
            glfwSetWindowShouldClose(window, GL_TRUE);

        switch(action) {
            case GLFW_PRESS:
                inputHandler.key(scancode, true);
                break;
            case GLFW_RELEASE:
                inputHandler.key(scancode, false);
                break;
            default:break;
        }
    });
    glfwSetCursorPosCallback(_window, [](GLFWwindow* window, double xpos, double ypos) {
        inputHandler.pos(xpos, ypos);
    });
    glfwSetMouseButtonCallback(_window, [](GLFWwindow* window, int button, int action, int mods) {
        switch(button) {
            case GLFW_MOUSE_BUTTON_LEFT:
                switch(action) {
                    case GLFW_PRESS:
                        inputHandler.leftDown(true);
                        break;
                    case GLFW_RELEASE:
                        inputHandler.leftDown(false);
                        break;
                    default:break;
                }
                break;
            case GLFW_MOUSE_BUTTON_MIDDLE:
                switch(action) {
                    case GLFW_PRESS:
                        inputHandler.wheelDown(true);
                        break;
                    case GLFW_RELEASE:
                        inputHandler.wheelDown(false);
                        break;
                    default:break;
                }
                break;
            case GLFW_MOUSE_BUTTON_RIGHT:
                switch(action) {
                    case GLFW_PRESS:
                        inputHandler.rightDown(true);
                        break;
                    case GLFW_RELEASE:
                        inputHandler.rightDown(false);
                        break;
                    default:break;
                }
                break;
            default:break;
        }
    });
    glfwSetScrollCallback(_window, [&](GLFWwindow* window, double xoffset, double yoffset) {
        inputHandler.delWheel(yoffset);
    });

    inputHandler.registerMouseListener([](InputHandler::MouseState &mouseState) {

    });
}

void Window::initloop() {
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glClearColor(0.2f, 0.2f, 0.2f, 1.f);

    while (!glfwWindowShouldClose(_window)) {
        int width, height;
        glfwGetFramebufferSize(_window, &width, &height);
        glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glBindVertexArray(vao);

        for (Painter* painter : _painters) {
            painter->draw();
        }

        glfwSwapBuffers(_window);
        glfwPollEvents();
    }

    glfwDestroyWindow(_window);
    glfwTerminate();
    exit(EXIT_SUCCESS);
}

void Window::handleMouseInput(InputHandler::MouseState &mouseState) {

}

void Window::addPainter(Painter *painter) {
    _painters.push_back(painter);
}

void Window::removePainter(Painter *painter) {
    _painters.erase(std::remove(_painters.begin(), _painters.end(), painter), _painters.end());
}
