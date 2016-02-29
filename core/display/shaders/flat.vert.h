//
// Created by austin on 2/29/16.
//

#ifndef FLUIDSOLVER_FLAT_VERT_H_H
#define FLUIDSOLVER_FLAT_VERT_H_H

const char* flat_vert = R"(
#version 150

in vec3 v_pos;
in vec3 v_col;

out vec3 f_col;

void main() {
    f_col = v_col;
    gl_Position = vec4(v_pos, 1);
}
)";

#endif //FLUIDSOLVER_FLAT_VERT_H_H