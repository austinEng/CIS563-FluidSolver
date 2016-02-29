//
// Created by austin on 2/29/16.
//

const char* flat_frag = R"(
#version 150

in vec3 f_col;
out vec4 out_Col;

void main() {
    out_Col = vec4(f_col.rgb, 1);
}
)";
