//
// Created by austin on 2/28/16.
//

#ifndef FLUIDSOLVER_GRIDATTR_FRAG_H
#define FLUIDSOLVER_GRIDATTR_FRAG_H

const char* gridAttr_frag = R"(
#version 150

uniform vec3 u_col;
out vec4 out_Col;
in float f_scale;

void main() {
    out_Col = vec4(f_scale*u_col.rgb, 1);
}
)";
#endif