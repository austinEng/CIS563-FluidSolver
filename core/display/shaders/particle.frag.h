//
// Created by austin on 2/28/16.
//

#ifndef FLUIDSOLVER_PARTICLE_FRAG_H
#define FLUIDSOLVER_PARTICLE_FRAG_H

const char* particle_frag = R"(
#version 150

in vec3 f_col;
out vec4 out_Col;

void main() {
    out_Col = vec4(f_col.rgb, 1);
}
)";
#endif