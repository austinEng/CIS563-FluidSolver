//
// Created by austin on 3/28/16.
//

#include "GridScalarAttributePainter.h"
#include <core/display/shaders/gridScal.frag.h>
#include <core/display/shaders/gridScal.vert.h>
#include <iostream>

GridScalarAttributePainter::GridScalarAttributePainter(
        Grid<float> *grid, float rangeStart, float rangeEnd,
        float ptSizeStart, float ptSizeEnd, const glm::vec3 &colorStart,
        const glm::vec3 &colorEnd)
        :
        GridScalarAttributePainter(grid, rangeStart, rangeEnd, ptSizeStart, ptSizeEnd, colorStart, colorEnd, FLOAT) {

}

GridScalarAttributePainter::GridScalarAttributePainter(
        Grid<int> *grid, float rangeStart, float rangeEnd,
        float ptSizeStart, float ptSizeEnd, const glm::vec3 &colorStart,
        const glm::vec3 &colorEnd)
        :
        GridScalarAttributePainter(grid, rangeStart, rangeEnd, ptSizeStart, ptSizeEnd, colorStart, colorEnd, INT) {

}

template <typename T>
GridScalarAttributePainter::GridScalarAttributePainter(Grid<T> *grid, float rangeStart, float rangeEnd,
                                                       float ptSizeStart, float ptSizeEnd, const glm::vec3 &colorStart,
                                                       const glm::vec3 &colorEnd, Type type) : type(type) {

    if (type == FLOAT) {
        _attributesF = &dynamic_cast<Grid<float>*>(grid)->_contents;
        MAX_ATTRIBUTES = _attributesF->size();
    } else if (type == INT) {
        _attributesI = &dynamic_cast<Grid<int>*>(grid)->_contents;
        MAX_ATTRIBUTES = _attributesI->size();
    }

    std::vector<int> indices;
    for (int i = 0; i < MAX_ATTRIBUTES; i++) {
        indices.push_back(i);
    }

    // compile shaders
    GLuint gridScalVert = compileShader(gridScal_vert, GL_VERTEX_SHADER);
    GLuint gridScalFrag = compileShader(gridScal_frag, GL_FRAGMENT_SHADER);

    std::vector<GLuint> programs = {gridScalVert, gridScalFrag};
    prog = makeProgram(programs);

    // setup shader locations
    unifViewProj = glGetUniformLocation(prog, "u_viewProj");
    if (type == FLOAT) {
        attrData = glGetAttribLocation(prog, "f_data");
    } else if (type == INT) {
        attrData = glGetAttribLocation(prog, "i_data");
    }

    unifColStart = glGetUniformLocation(prog, "u_colStart");
    unifColEnd = glGetUniformLocation(prog, "u_colEnd");
    unifSizeStart = glGetUniformLocation(prog, "u_sizeStart");
    unifSizeEnd = glGetUniformLocation(prog, "u_sizeEnd");
    unifRangeStart = glGetUniformLocation(prog, "u_rangeStart");
    unifRangeEnd = glGetUniformLocation(prog, "u_rangeEnd");
    unifType = glGetUniformLocation(prog, "u_type");

    unifCellSize = glGetUniformLocation(prog, "u_cellSize");
    unifCellCount = glGetUniformLocation(prog, "u_cellCount");
    unifOrigin = glGetUniformLocation(prog, "u_origin");

    // make a buffer for the indices
    glGenBuffers(1, &index_buffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, MAX_ATTRIBUTES * sizeof(int), indices.data(), GL_STATIC_DRAW);

    // make a buffer for the attributes
    if (type == FLOAT) {
        glGenBuffers(1, &attribute_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, attribute_buffer);
        glBufferData(GL_ARRAY_BUFFER, MAX_ATTRIBUTES * sizeof(float), NULL, GL_STREAM_DRAW);

        glUniform1i(unifType, 0);
    } else if (type == INT) {
        glGenBuffers(1, &attribute_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, attribute_buffer);
        glBufferData(GL_ARRAY_BUFFER, MAX_ATTRIBUTES * sizeof(int), NULL, GL_STREAM_DRAW);

        glUniform1i(unifType, 1);
    }

    // set grid uniforms
    glUseProgram(prog);
    glUniform1f(unifCellSize, grid->_cellSize);
    glm::ivec3 count = glm::ivec3(grid->_countX, grid->_countY, grid->_countZ);
    glUniform3iv(unifCellCount, 1, &(count.x));
    glm::vec3 o = grid->_origin + grid->_offset;
    glUniform3fv(unifOrigin, 1, &(o[0]));

    glUniform3fv(unifColStart, 1, &(colorStart[0]));
    glUniform3fv(unifColEnd, 1, &(colorEnd[0]));
    glUniform1f(unifRangeStart, rangeStart);
    glUniform1f(unifRangeEnd, rangeEnd);
    glUniform1f(unifSizeStart, ptSizeStart);
    glUniform1f(unifSizeEnd, ptSizeEnd);
}

void GridScalarAttributePainter::draw() const {
    if (type == FLOAT) {
        if (_attributesF != nullptr) {
//            for (float f : *_attributesF) {
//                std:: cout << f << std::endl;
//            }

            glUseProgram(prog);
            glEnable(GL_PROGRAM_POINT_SIZE);

            // bind and send new data
            glBindBuffer(GL_ARRAY_BUFFER, attribute_buffer);
            glBufferData(GL_ARRAY_BUFFER, MAX_ATTRIBUTES * sizeof(float), NULL, GL_STREAM_DRAW);
            glBufferSubData(GL_ARRAY_BUFFER, 0, MAX_ATTRIBUTES * sizeof(float), _attributesF->data());

            glEnableVertexAttribArray(attrData);
            glVertexAttribPointer(attrData, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
            glDrawElements(GL_POINTS, _attributesF->size(), GL_UNSIGNED_INT, 0);

            glDisableVertexAttribArray(attrData);
        }
    } else if (type == INT) {
        if (_attributesI != nullptr) {
            glUseProgram(prog);
            glEnable(GL_PROGRAM_POINT_SIZE);

            // bind and send new data
            glBindBuffer(GL_ARRAY_BUFFER, attribute_buffer);
            glBufferData(GL_ARRAY_BUFFER, MAX_ATTRIBUTES * sizeof(int), NULL, GL_STREAM_DRAW);
            glBufferSubData(GL_ARRAY_BUFFER, 0, MAX_ATTRIBUTES * sizeof(int), _attributesI->data());

            glEnableVertexAttribArray(attrData);
            glVertexAttribIPointer(attrData, 1, GL_INT, sizeof(int), (void*)0);

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
            glDrawElements(GL_POINTS, _attributesI->size(), GL_UNSIGNED_INT, 0);

            glDisableVertexAttribArray(attrData);
        }
    }
}

void GridScalarAttributePainter::setViewProj(const float *viewProj) {
    glUseProgram(prog);
    glUniformMatrix4fv(unifViewProj, 1, GL_FALSE, viewProj);
}
