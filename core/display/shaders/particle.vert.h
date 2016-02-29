//
// Created by austin on 2/28/16.
//

const char* particle_vert = R"(
#version 150

in vec3 v_pos;
in vec3 v_vel;
in vec3 v_col;

out vec3 f_col;

void main() {
    f_col = v_col;
    gl_Position = vec4(v_pos, 1);
}
)";