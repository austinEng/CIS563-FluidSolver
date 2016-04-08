//
// Created by austin on 3/28/16.
//

#ifndef FLUIDSOLVER_GRIDSCALARATTRIBUTEPAINTER_H
#define FLUIDSOLVER_GRIDSCALARATTRIBUTEPAINTER_H

#include "Painter.h"
#include <core/solver/grid/Grid.h>

class GridScalarAttributePainter : public Painter {
public:
    enum Type {
        INT,
        FLOAT
    };

    template <typename T> explicit GridScalarAttributePainter(Grid<T>* grid,
                                        float rangeStart,
                                        float rangeEnd,
                                        float ptSizeStart,
                                        float ptSizeEnd,
                                        const glm::vec3 &colorStart,
                                        const glm::vec3 &colorEnd,
                                        Type type);

    explicit GridScalarAttributePainter(Grid<float>* grid,
                                        float rangeStart,
                                        float rangeEnd,
                                        float ptSizeStart,
                                        float ptSizeEnd,
                                        const glm::vec3 &colorStart,
                                        const glm::vec3 &colorEnd);

    explicit GridScalarAttributePainter(Grid<int>* grid,
                                        float rangeStart,
                                        float rangeEnd,
                                        float ptSizeStart,
                                        float ptSizeEnd,
                                        const glm::vec3 &colorStart,
                                        const glm::vec3 &colorEnd);

    virtual void draw() const;
    virtual void setViewProj(const float* viewProj);

private:
    Type type;
    unsigned int MAX_ATTRIBUTES = 10000;
    GLuint index_buffer;
    GLuint attribute_buffer;

    GLint attrIndex;
    GLint attrData;

    GLint unifViewProj;

    GLint unifColStart;
    GLint unifColEnd;
    GLint unifSizeStart;
    GLint unifSizeEnd;
    GLint unifRangeStart;
    GLint unifRangeEnd;
    GLint unifType;

    GLint unifCellSize;
    GLint unifCellCount;
    GLint unifOrigin;

    std::vector<float>* _attributesF;
    std::vector<int>* _attributesI;
};


#endif //FLUIDSOLVER_GRIDSCALARATTRIBUTEPAINTER_H
