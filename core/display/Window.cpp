//
// Created by austin on 2/25/16.
//

#include "Window.h"
#include <iostream>

#include "shaders/particle.frag.h"
#include "shaders/particle.vert.h"

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

GLuint compileShader(const char* shader, GLenum type) {
    GLuint shaderId = glCreateShader(type);
    glShaderSourceARB(shaderId, 1, &shader, NULL);
    glCompileShader(shaderId);

    GLint success = 0;
    glGetShaderiv(shaderId, GL_COMPILE_STATUS, &success);
    if (success == GL_FALSE) {
        fprintf(stderr, "Failed to compile shader!\n%s\n", shader);

        GLint maxLength = 0;
        glGetShaderiv(shaderId, GL_INFO_LOG_LENGTH, &maxLength);

        std::vector<GLchar> errorLog(maxLength);
        glGetShaderInfoLog(shaderId, maxLength, &maxLength, &errorLog[0]);

        fprintf(stderr, "%s\n", &errorLog[0]);

        glDeleteShader(shaderId);
    }

    return shaderId;
}

GLuint makeProgram(std::vector<GLuint> programs) {
    GLuint prog = glCreateProgram();
    for (GLuint program : programs) {
        glAttachShader(prog, program);
    }
    glLinkProgram(prog);

    for (GLuint program : programs) {
        glDetachShader(prog, program);
        glDeleteShader(program);
    }

    GLint linked;
    glGetProgramiv(prog, GL_LINK_STATUS, &linked);
    if (!linked) {
        std::cerr << "Failed to link program!" << std::endl;

        GLint length;
        glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &length);
        if ( length > 0 ){
            std::vector<char> ProgramErrorMessage(length+1);
            glGetProgramInfoLog(prog, length, NULL, &ProgramErrorMessage[0]);
            fprintf(stderr, "%s\n", &ProgramErrorMessage[0]);
        }
    }
    return prog;
}

void Window::initloop() {
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    GLuint particleVert = compileShader(particle_vert.c_str(), GL_VERTEX_SHADER);
    GLuint particleFrag = compileShader(particle_frag.c_str(), GL_FRAGMENT_SHADER);

    std::vector<GLuint> programs = {particleVert, particleFrag};
    GLuint particleProg = makeProgram(programs);

    GLint attrPos = glGetAttribLocation(particleProg, "v_pos");;
    GLint attrVel = glGetAttribLocation(particleProg, "v_vel");
    GLint attrCol = glGetAttribLocation(particleProg, "v_col");

    glClearColor(0.2f, 0.2f, 0.2f, 1.f);

    int MaxParticles = 1000;

    GLuint particles_position_buffer;
    glGenBuffers(1, &particles_position_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
    //glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);
    glBufferData(GL_ARRAY_BUFFER, MaxParticles * sizeof(glm::vec3), NULL, GL_STREAM_DRAW);

    GLuint particles_color_buffer;
    glGenBuffers(1, &particles_color_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
    //glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW);
    glBufferData(GL_ARRAY_BUFFER, MaxParticles * sizeof(glm::vec3), NULL, GL_STREAM_DRAW);

    GLuint particles_buffer;
    glGenBuffers(1, &particles_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particles_buffer);
    //glBufferData(GL_ARRAY_BUFFER, MaxParticles * sizeof(Fluid::Particle), NULL, GL_STREAM_DRAW);
    glBufferData(GL_ARRAY_BUFFER, MaxParticles * sizeof(glm::vec3), NULL, GL_STREAM_DRAW);


    std::vector<Fluid::Particle> particles;
    Fluid::Particle p1;
    p1.pos = glm::vec3(0, 0, 0);

    Fluid::Particle p2;
    p2.pos = glm::vec3(0, 0.8, 0);

    Fluid::Particle p3;
    p3.pos = glm::vec3(0.8, 0.8, 0);

    particles.push_back(p1);
    particles.push_back(p2);
    particles.push_back(p3);

    while (!glfwWindowShouldClose(_window)) {
        int width, height;
        glfwGetFramebufferSize(_window, &width, &height);
        glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glBindBuffer(GL_ARRAY_BUFFER, particles_buffer);
        glBufferData(GL_ARRAY_BUFFER, MaxParticles * sizeof(Fluid::Particle), NULL, GL_STREAM_DRAW);
        glBufferSubData(GL_ARRAY_BUFFER, 0, MaxParticles * sizeof(Fluid::Particle), &particles[0]);

        glEnableVertexAttribArray(attrPos);
        glVertexAttribPointer(attrPos, 3, GL_FLOAT, GL_FALSE, sizeof(Fluid::Particle), (void*)offsetof(Fluid::Particle, pos));

        glEnableVertexAttribArray(attrVel);
        glVertexAttribPointer(attrVel, 3, GL_FLOAT, GL_FALSE, sizeof(Fluid::Particle), (void*)offsetof(Fluid::Particle, vel));

        glEnableVertexAttribArray(attrCol);
        glVertexAttribPointer(attrCol, 3, GL_FLOAT, GL_FALSE, sizeof(Fluid::Particle), (void*)offsetof(Fluid::Particle, col));

        glUseProgram(particleProg);
        glBindVertexArray(vao);

        glPointSize(3);

        glDrawArrays(GL_POINTS, 0, (GLsizei) particles.size());

        glDisableVertexAttribArray(attrPos);
        glDisableVertexAttribArray(attrVel);
        glDisableVertexAttribArray(attrCol);

        glfwSwapBuffers(_window);
        glfwPollEvents();
    }

    glfwDestroyWindow(_window);
    glfwTerminate();
    exit(EXIT_SUCCESS);
}

void Window::handleMouseInput(InputHandler::MouseState &mouseState) {

}
