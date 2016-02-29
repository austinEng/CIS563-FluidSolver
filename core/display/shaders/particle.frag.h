//
// Created by austin on 2/28/16.
//

std::string particle_frag = R"(
#version 150

in vec3 f_col;
out vec4 out_Col;

void main() {
    out_Col = vec4(f_col.rgb, 1);
    if (f_col.b == 0) {
        //out_Col = vec4(1,1,1,1);
    }
}
)";