//
// Created by austin on 3/28/16.
//

#ifndef FLUIDSOLVER_GRIDSCAL_FRAG_H
#define FLUIDSOLVER_GRIDSCAL_FRAG_H

const char* gridScal_frag = R"(
#version 150

uniform vec3 u_colStart;
uniform vec3 u_colEnd;

in float amount;
out vec4 out_Col;

void main() {
    out_Col = vec4((1 - amount)*u_colStart + amount*u_colEnd, 1);
}
)";

#endif //FLUIDSOLVER_GRIDSCAL_FRAG_H
