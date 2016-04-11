//
// Created by austin on 2/28/16.
//

#ifndef FLUIDSOLVER_PARTICLE_FRAG_H
#define FLUIDSOLVER_PARTICLE_FRAG_H

const char* particle_frag = R"(
#version 150

in vec3 f_col;
in vec3 f_vel;
out vec4 out_Col;

void main() {
    float fac = length(f_vel) / 20.0;
    vec3 col = f_col + vec3(fac,fac,fac);
    out_Col = vec4(col, 1);
}
)";
#endif