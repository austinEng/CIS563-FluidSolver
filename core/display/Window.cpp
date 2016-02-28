//
// Created by austin on 2/25/16.
//

#include "Window.h"
#include <stdio.h>
#include <stdlib.h>
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
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); //We don't want the old OpenGL

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

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
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

    glClearColor(0.2f, 0.2f, 0.2f, 1.f);

    static const GLfloat particle_vertex_data[] = {
        -0.5f, -0.5f, 0.0f,
        0.5f, -0.5f, 0.0f,
        -0.5f, 0.5f, 0.0f,
        0.5f, 0.5f, 0.0f,
    };

    GLuint particles_position_buffer;

    glGenBuffers(1, &particles_position_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);

    // Initialize with empty (NULL) buffer : it will be updated later, each frame.
    glBufferData(GL_ARRAY_BUFFER, 1000 * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);

    while (!glfwWindowShouldClose(_window)) {
        float ratio;
        int width, height;
        glfwGetFramebufferSize(_window, &width, &height);
        ratio = width / (float) height;

        glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-ratio, ratio, -1.f, 1.f, 1.f, -1.f);
        glMatrixMode(GL_MODELVIEW);

        glLoadIdentity();
        glRotatef((float) glfwGetTime() * 50.f, 0.f, 0.f, 1.f);
        glBegin(GL_TRIANGLES);
        glColor3f(1.f, 0.f, 0.f);
        glVertex3f(-0.6f, -0.4f, 0.f);
        glColor3f(0.f, 1.f, 0.f);
        glVertex3f(0.6f, -0.4f, 0.f);
        glColor3f(0.f, 0.f, 1.f);
        glVertex3f(0.f, 0.6f, 0.f);
        glEnd();

        glfwSwapBuffers(_window);
        glfwPollEvents();
    }

    glfwDestroyWindow(_window);
    glfwTerminate();
    exit(EXIT_SUCCESS);
}

void Window::handleMouseInput(InputHandler::MouseState &mouseState) {

}
