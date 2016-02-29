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
Window::Window(int w, int h, const char* title) : _window(nullptr), camera(w, h), _w(w), _h(h) {
    glfwSetErrorCallback(error_callback);

    if (!glfwInit()) exit(EXIT_FAILURE);

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

    _window = glfwCreateWindow(w, h, title, NULL, NULL);

    if (!_window) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(_window);
    glfwSwapInterval(1);
    setupInputCBs();

    glewExperimental= GL_TRUE;
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        exit(EXIT_FAILURE);
    }
}

Window::~Window() {

}

void Window::setupInputCBs() {

    glfwSetKeyCallback(_window, [](GLFWwindow *window, int key, int scancode, int action, int mods) {
        switch(action) {
            case GLFW_PRESS:
                inputHandler.key(key, true);
                break;
            case GLFW_RELEASE:
                inputHandler.key(key, false);
                break;
            default: break;
        }
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

    glfwSetWindowSizeCallback(_window, [](GLFWwindow *window, int width, int height) {
        inputHandler.windowResized(width, height);
    });

    inputHandler.registerMouseListener([&](InputHandler::MouseState &mouseState) {
        if (!fequal(mouseState.delWheel, 0.0)) {
            glm::vec3 vec = camera.tgt - camera.eye;
            float fac = glm::min(glm::abs(glm::length(vec)/5.f), 1.f);
            camera.zoom -= (float)mouseState.delWheel * fac;
            camera.recomputeEye();
            updateCamera();
        }

        if (mouseState.wheelDragging) {
            if (inputHandler.key(340)) {
                float x = (float) (_w/2 - mouseState.delX);
                float y = (float) (_h/2 - mouseState.delY);

                float sx = (2*x / _w) - 1.f;
                float sy = 1.f - (2*y / _h);

                float alpha = camera.fovy / 2;
                float len = glm::length(camera.tgt - camera.eye);
                glm::vec3 V = camera.up*(float)(len*tan(alpha));
                glm::vec3 H = camera.right*(float)(len*(_w / _h)*tan(alpha));

                camera.tgt = camera.tgt + sx*H + sy*V;
                camera.recomputeEye();
                updateCamera();
                return;
            }
            glm::vec4 y = glm::vec4(0,1,0,0);
            glm::vec4 diff(mouseState.delX / _w, mouseState.delY / _h, 0, 0);
            float a = (float) acos(glm::dot(y, diff) / (glm::length(y) * (glm::length(diff))));
            glm::vec4 para;
            if (diff[0] > 0) {
                para = glm::mat4_cast(glm::angleAxis(-a, glm::vec3(camera.look[0], camera.look[1], camera.look[2]))) * glm::vec4(camera.up, 1);
            } else {
                para = glm::mat4_cast(glm::angleAxis(a, glm::vec3(camera.look[0], camera.look[1], camera.look[2]))) * glm::vec4(camera.up, 1);
            }
            glm::vec3 perp = glm::normalize(glm::cross(camera.look, glm::vec3(para)));

            glm::mat4 rot = glm::mat4_cast(glm::angleAxis(-2*PI*glm::length(diff), perp));
            camera.rotation = rot * camera.rotation;
            camera.up = glm::vec3(rot * glm::vec4(camera.up, 0));
            camera.right = glm::vec3(rot * glm::vec4(camera.right, 0));
            camera.look = glm::vec3(rot * glm::vec4(camera.look, 0));
            camera.recomputeEye();
            updateCamera();
        }
    });

    inputHandler.registerWindowListener([&](int w, int h){
        _w = w;
        _h = h;
        glViewport(0, 0, w, h);
        camera.resize(w, h);
        updateCamera();
    });
}

void Window::initloop(std::function<void(void)> predraw) {
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    for (Painter* painter : _painters) {
        painter->setViewProj(glm::value_ptr(camera.viewProj()));
    }

    glClearColor(0.2f, 0.2f, 0.2f, 1.f);
    glEnable(GL_DEPTH_TEST);

    while (!glfwWindowShouldClose(_window)) {
        predraw();
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

void Window::updateCamera() {
    camera.recompute();
    for (Painter* painter : _painters) {
        painter->setViewProj(glm::value_ptr(camera.viewProj()));
    }
}
