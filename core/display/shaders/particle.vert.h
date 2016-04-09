//
// Created by austin on 2/28/16.
//

#ifndef FLUIDSOLVER_PARTICLE_VERT_H
#define FLUIDSOLVER_PARTICLE_VERT_H

const char* particle_vert = R"(
#version 150

uniform mat4 u_viewProj;

in vec3 v_pos;
in vec3 v_vel;
in vec3 v_col;

out vec3 f_col;
out vec3 f_vel;

void main() {
    f_col = v_col;
    f_vel = v_vel;
    gl_Position = u_viewProj * vec4(v_pos, 1);
    gl_PointSize = 3;
}
)";

#endif