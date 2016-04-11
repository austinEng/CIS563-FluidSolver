//
// Created by austin on 2/25/16.
//

#include "Window.h"
#include <algorithm>
#include <iostream>
#include <fstream>

static void error_callback(int error, const char* description) {
    fputs(description, stderr);
}

InputHandler &inputHandler = InputHandler::getInputHandler();

Window::Window(const char *title) : Window(1200, 800, title) { }
Window::Window(int w, int h) : Window(w, h, "GL Window"){ }
Window::Window(int w, int h, const char* title) : _window(nullptr), camera(w, h), _w(w), _h(h),
    loadSceneCB(NULL) {
    glfwSetErrorCallback(error_callback);

    if (!glfwInit()) exit(EXIT_FAILURE);

    // use antialiasing
    glfwWindowHint(GLFW_SAMPLES, 4);

    // set version to OpenGL 3.3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    #ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

    _window = glfwCreateWindow(w, h, title, NULL, NULL);
    pixels.resize(_w*_h*4);

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

    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);

    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, texture, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

Window::~Window() {

}

void Window::initializeTweakBar() {
    int w, h;
    glfwGetWindowSize(_window, &w, &h);
    TwInit(TW_OPENGL_CORE, NULL);
    TwWindowSize(w, h);
    TwBar *myBar;
    myBar = TwNewBar("Settings");
    TwAddButton(myBar, "loadsceneBtn", loadSceneCB, NULL, " label='Load Scene'");
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

        TwEventKeyGLFW(key, action);
        TwEventCharGLFW(key, action);
    });
    glfwSetCursorPosCallback(_window, [](GLFWwindow* window, double xpos, double ypos) {
        inputHandler.pos(xpos, ypos);
        TwEventMousePosGLFW(xpos, ypos);
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
        TwEventMouseButtonGLFW(button, action);
    });
    glfwSetScrollCallback(_window, [&](GLFWwindow* window, double xoffset, double yoffset) {
        inputHandler.delWheel(yoffset);
        TwEventMouseWheelGLFW(yoffset);
    });

    glfwSetWindowSizeCallback(_window, [](GLFWwindow *window, int width, int height) {
        inputHandler.windowResized(width, height);
        TwWindowSize(width, height);
    });

    inputHandler.registerMouseListener([&](InputHandler::MouseState &mouseState) {
        if (!fequal(mouseState.delWheel, 0.0)) {
            // change camera zoom level based on scroll direction
            glm::vec3 vec = camera.tgt - camera.eye;
            // limit zoom when very near to target
            float fac = glm::min(glm::abs(glm::length(vec)/5.f), 1.f);
            camera.zoom -= (float)mouseState.delWheel * fac;
            camera.recomputeEye();
            updateCamera();
        }

        if (mouseState.wheelDragging) {
            if (inputHandler.key(340)) {
                // pixel position offset from center
                float x = (float) (_w/2 - mouseState.delX);
                float y = (float) (_h/2 - mouseState.delY);

                // offset in ndc
                float sx = (2*x / _w) - 1.f;
                float sy = 1.f - (2*y / _h);

                // project camera up amd right axes
                float alpha = camera.fovy / 2;
                float len = glm::length(camera.tgt - camera.eye);
                glm::vec3 V = camera.up*(float)(len*tan(alpha));
                glm::vec3 H = camera.right*(float)(len*(_w / _h)*tan(alpha));

                camera.tgt = camera.tgt + sx*H + sy*V;
                camera.recomputeEye();
                updateCamera();
                return;
            }

            glm::vec4 y = glm::vec4(0,1,0,0); // y axis vector
            glm::vec4 diff(mouseState.delX / _w, mouseState.delY / _h, 0, 0); // mouse offset
            float a = (float) acos(glm::dot(y, diff) / (glm::length(y) * (glm::length(diff)))); // calculate offset angle from y axis
            glm::vec4 para; // parallel axis to mouse movement
            if (diff[0] > 0) {
                para = glm::mat4_cast(glm::angleAxis(-a, glm::vec3(camera.look[0], camera.look[1], camera.look[2]))) * glm::vec4(camera.up, 1);
            } else {
                para = glm::mat4_cast(glm::angleAxis(a, glm::vec3(camera.look[0], camera.look[1], camera.look[2]))) * glm::vec4(camera.up, 1);
            }
            glm::vec3 perp = glm::normalize(glm::cross(camera.look, glm::vec3(para))); // perpendicular axis to mouse movement

            // rotate camera on perpendicular axis
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

//https://danielbeard.wordpress.com/2011/06/06/image-saving-code-c/
void Window::saveImage(const std::string &filename) {
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glViewport(0,0,_w,_h);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    for (Painter* painter : _painters) {
        painter->draw();
    }
//    std::vector<char> pixels(this->_w*this->_h*4);

    glReadPixels(0,0,_w,_h, GL_BGRA, GL_UNSIGNED_BYTE, &(pixels[0]));

    std::ofstream o(filename.c_str(), std::ios::out | std::ios::binary);
    o.put(0);
    o.put(0);
    o.put(2);                         /* uncompressed RGB */
    o.put(0); 		o.put(0);
    o.put(0); 	o.put(0);
    o.put(0);
    o.put(0); 	o.put(0);           /* X origin */
    o.put(0); 	o.put(0);           /* y origin */
    o.put((_w & 0x00FF));
    o.put((_w & 0xFF00) / 256);
    o.put((_h & 0x00FF));
    o.put((_h & 0xFF00) / 256);
    o.put(32);                        /* 24 bit bitmap */
    o.put(0);

    for (int i=0;i<_w*_h*4;i+=4) {
//        std::cout << (unsigned int)pixels[i] << "," << (unsigned int)pixels[i+1] << "," << (unsigned int)pixels[i+2] << "," << (unsigned int)pixels[i+3] << std::endl;
        o.put(pixels[i+0]);
        o.put(pixels[i+1]);
        o.put(pixels[i+2]);
        o.put(255);
    }

    o.close();

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
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


//        TwDraw();

        glfwSwapBuffers(_window);
        glfwPollEvents();
    }

    TwTerminate();

    glfwDestroyWindow(_window);
    glfwTerminate();
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
    // send camera uniforms to painters
    for (Painter* painter : _painters) {
        painter->setViewProj(glm::value_ptr(camera.viewProj()));
    }
}
