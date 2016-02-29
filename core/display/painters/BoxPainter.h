//
// Created by austin on 2/29/16.
//

#ifndef FLUIDSOLVER_BOXPAINTER_H
#define FLUIDSOLVER_BOXPAINTER_H

#include "Painter.h"
#include <core/geometry/Box.h>

class BoxPainter : public Painter {
public:
    BoxPainter(Box* box);
    void update();
    void draw() const;

private:
    GLuint vertex_buffer;
    GLuint index_buffer;

    GLint attrPos;
    GLint attrCol;

    Box* _box;

    void create();
    void destroy();
};


#endif //FLUIDSOLVER_BOXPAINTER_H
