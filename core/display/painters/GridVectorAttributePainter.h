//
// Created by austin on 3/22/16.
//

#ifndef FLUIDSOLVER_GRIDVECTORATTRIBUTEPAINTER_H
#define FLUIDSOLVER_GRIDVECTORATTRIBUTEPAINTER_H

#include "Painter.h"
#include <core/solver/grid/Grid.h>

class GridVectorAttributePainter : public Painter {
public:
    explicit GridVectorAttributePainter(Grid<float>* grid, float ptSize, const glm::vec3 &color, const glm::vec3 &dir);
    virtual void draw() const;
    virtual void setViewProj(const float* viewProj);

private:
    unsigned int MAX_ATTRIBUTES = 10000;
    GLuint index_buffer;
    GLuint attribute_buffer;

    GLint unifViewProj;
    GLint attrIndex;
    GLint attrData;
    GLint unifCol;
    GLint unifCellSize;
    GLint unifCellCount;
    GLint unifOrigin;
    GLint unifVec;

    GLfloat _ptSize;
    std::vector<float>* _attributes;
};


#endif //FLUIDSOLVER_GRIDATTRIBUTEPAINTER_H
