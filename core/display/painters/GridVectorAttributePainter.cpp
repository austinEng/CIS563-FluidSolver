//
// Created by austin on 3/22/16.
//

#include "GridVectorAttributePainter.h"
#include <core/display/shaders/gridAttr.frag.h>
#include <core/display/shaders/gridAttr.vert.h>
#include <core/display/shaders/gridAttr.geo.h>

GridVectorAttributePainter::GridVectorAttributePainter(Grid<float> *grid,
                                                       float ptSize,
                                                       const glm::vec3 &color,
                                                       const glm::vec3 &dir) :
        _ptSize(ptSize), _attributes(&grid->_contents) {
    MAX_ATTRIBUTES = (unsigned int) _attributes->size();

    std::vector<int> indices;
    for (int i = 0; i < _attributes->size(); i++) {
        indices.push_back(i);
    }

    // compile shaders
    GLuint gridAttrVert = compileShader(gridAttr_vert, GL_VERTEX_SHADER);
    GLuint gridAttrGeo = compileShader(gridAttr_geo, GL_GEOMETRY_SHADER);
    GLuint gridAttrFrag = compileShader(gridAttr_frag, GL_FRAGMENT_SHADER);

    std::vector<GLuint> programs = {gridAttrVert, gridAttrGeo, gridAttrFrag};
    prog = makeProgram(programs);

    // setup shader locations
    unifViewProj = glGetUniformLocation(prog, "u_viewProj");
    attrData = glGetAttribLocation(prog, "v_data");
    unifCol = glGetUniformLocation(prog, "u_col");
    unifCellSize = glGetUniformLocation(prog, "u_cellSize");
    unifCellCount = glGetUniformLocation(prog, "u_cellCount");
    unifOrigin = glGetUniformLocation(prog, "u_origin");
    unifVec = glGetUniformLocation(prog, "u_vec");

    // make a buffer for the indices
    glGenBuffers(1, &index_buffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, MAX_ATTRIBUTES * sizeof(int), indices.data(), GL_STATIC_DRAW);

    // make a buffer for the attributes
    glGenBuffers(1, &attribute_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, attribute_buffer);
    glBufferData(GL_ARRAY_BUFFER, MAX_ATTRIBUTES * sizeof(float), NULL, GL_STREAM_DRAW);

    // set grid uniforms
    glUseProgram(prog);
    glUniform1f(unifCellSize, grid->_cellSize);
    glm::ivec3 count = glm::ivec3(grid->_countX, grid->_countY, grid->_countZ);
    glUniform3iv(unifCellCount, 1, &(count.x));
    glm::vec3 o = grid->_origin + grid->_offset;
    glUniform3fv(unifOrigin, 1, &(o[0]));
    glUniform3fv(unifCol, 1, &(color[0]));
    glUniform3fv(unifVec, 1, &(dir[0]));
}

void GridVectorAttributePainter::draw() const {
    if (_attributes != nullptr) {
        glUseProgram(prog);

        // bind and send new data
        glBindBuffer(GL_ARRAY_BUFFER, attribute_buffer);
        glBufferData(GL_ARRAY_BUFFER, MAX_ATTRIBUTES * sizeof(float), NULL, GL_STREAM_DRAW);
        glBufferSubData(GL_ARRAY_BUFFER, 0, MAX_ATTRIBUTES * sizeof(float), _attributes->data());

        glEnableVertexAttribArray(attrData);
        glVertexAttribPointer(attrData, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);

        glPointSize(_ptSize);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
        glDrawElements(GL_POINTS, _attributes->size(), GL_UNSIGNED_INT, 0);

        glDisableVertexAttribArray(attrData);
    }
}

void GridVectorAttributePainter::setViewProj(const float *viewProj) {
    glUseProgram(prog);
    glUniformMatrix4fv(unifViewProj, 1, GL_FALSE, viewProj);
}
