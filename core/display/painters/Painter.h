//
// Created by austin on 2/29/16.
//

#ifndef FLUIDSOLVER_PAINTER_H
#define FLUIDSOLVER_PAINTER_H

#include <GL/glew.h>
#include <vector>

class Painter {
public:
    virtual void draw() const = 0;
    virtual void setViewProj(const float* viewProj);

protected:
    GLuint prog;

    GLuint compileShader(const char* shader, GLenum type);
    GLuint makeProgram(const std::vector<GLuint> &programs);
};


#endif //FLUIDSOLVER_PAINTER_H
