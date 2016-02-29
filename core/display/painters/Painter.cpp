//
// Created by austin on 2/29/16.
//

#include "Painter.h"
#include <iostream>

GLuint Painter::compileShader(const char* shader, GLenum type) {
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

GLuint Painter::makeProgram(std::vector<GLuint> programs) {
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

void Painter::setViewProj(const float *viewProj) {

}
