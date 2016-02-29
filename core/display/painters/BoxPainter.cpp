//
// Created by austin on 2/29/16.
//

#include "BoxPainter.h"
#include <core/display/shaders/flat.vert.h>
#include <core/display/shaders/flat.frag.h>

struct vert {
    glm::vec3 pos;
    glm::vec3 col;
};

BoxPainter::BoxPainter(Box *box) : _box(box) {
    GLuint vert = compileShader(flat_vert, GL_VERTEX_SHADER);
    GLuint frag = compileShader(flat_frag, GL_FRAGMENT_SHADER);

    std::vector<GLuint> programs = {vert, frag};
    prog = makeProgram(programs);

    unifViewProj = glGetUniformLocation(prog, "u_viewProj");
    attrPos = glGetAttribLocation(prog, "v_pos");
    attrCol = glGetAttribLocation(prog, "v_col");

    glGenBuffers(1, &vertex_buffer);
    glGenBuffers(1, &index_buffer);

    create();
}

void BoxPainter::update() {
    if (_box != nullptr) {
        create();
        destroy();
    }
}

void BoxPainter::draw() const {
    if (_box != nullptr) {
        glUseProgram(prog);

        glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);

        glEnableVertexAttribArray(attrPos);
        glVertexAttribPointer(attrPos, 3, GL_FLOAT, GL_FALSE, sizeof(vert), (void*)offsetof(vert, pos));

        glEnableVertexAttribArray(attrCol);
        glVertexAttribPointer(attrCol, 3, GL_FLOAT, GL_FALSE, sizeof(vert), (void*)offsetof(vert, col));

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
        glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, (void*)0);

        glDisableVertexAttribArray(attrPos);
        glDisableVertexAttribArray(attrCol);
    }
}

void BoxPainter::create() {
    if (_box != nullptr) {
        vert verts[8];
        GLuint indices[24];
        
        float minX = _box->minX();
        float maxX = _box->maxX();
        float minY = _box->minY();
        float maxY = _box->maxY();
        float minZ = _box->minZ();
        float maxZ = _box->maxZ();
        
        verts[0].pos = glm::vec3(minX, minY, minZ);
        verts[1].pos = glm::vec3(minX, maxY, minZ);
        verts[2].pos = glm::vec3(minX, maxY, maxZ);
        verts[3].pos = glm::vec3(minX, minY, maxZ);
        verts[4].pos = glm::vec3(maxX, minY, minZ);
        verts[5].pos = glm::vec3(maxX, maxY, minZ);
        verts[6].pos = glm::vec3(maxX, maxY, maxZ);
        verts[7].pos = glm::vec3(maxX, minY, maxZ);

        indices[0] = 0; indices[1] = 1;
        indices[2] = 1; indices[3] = 2;
        indices[4] = 2; indices[5] = 3;
        indices[6] = 3; indices[7] = 0;

        indices[8] = 4; indices[9] = 5;
        indices[10] = 5; indices[11] = 6;
        indices[12] = 6; indices[13] = 7;
        indices[14] = 7; indices[15] = 4;

        indices[16] = 0; indices[17] = 4;
        indices[18] = 1; indices[19] = 5;
        indices[20] = 2; indices[21] = 6;
        indices[22] = 3; indices[23] = 7;

        for (int i = 0; i < 8; ++i) {
            verts[i].col = glm::vec3(1,1,1);
        }

        glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
        glBufferData(GL_ARRAY_BUFFER, 8 * sizeof(vert), verts, GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 24 * sizeof(GLuint), indices, GL_STATIC_DRAW);
    }
}

void BoxPainter::destroy() {
    glDeleteBuffers(1, &vertex_buffer);
    glDeleteBuffers(1, &index_buffer);
}

void BoxPainter::setViewProj(const float *viewProj) {
    glUseProgram(prog);
    glUniformMatrix4fv(unifViewProj, 1, GL_FALSE, viewProj);
}
