//
// Created by austin on 3/28/16.
//

#ifndef FLUIDSOLVER_GRIDSCAL_VERT_H
#define FLUIDSOLVER_GRIDSCAL_VERT_H

const char* gridScal_vert = R"(
#version 150

uniform mat4 u_viewProj;

uniform float u_cellSize;
uniform ivec3 u_cellCount;
uniform vec3 u_origin;

uniform float u_sizeStart;
uniform float u_sizeEnd;
uniform float u_rangeStart;
uniform float u_rangeEnd;

uniform int u_type;

in float f_data;
in int i_data;
out float amount;

void main() {

    int i = int(mod(gl_VertexID, u_cellCount.z));
    int j = int(mod(gl_VertexID / u_cellCount.z, u_cellCount.y));
    int k = int(gl_VertexID / (u_cellCount.y * u_cellCount.z));

    vec3 pos = vec3(float(i), float(j), float(k)) * u_cellSize + u_origin;

    gl_Position = u_viewProj * vec4(pos, 1);

    if (u_type > 0) {
        amount = (float(i_data) - u_rangeStart) / (u_rangeEnd - u_rangeStart);
    } else {
        amount = (f_data - u_rangeStart) / (u_rangeEnd - u_rangeStart);
    }
    amount = (float(i_data) - u_rangeStart) / (u_rangeEnd - u_rangeStart);

    gl_PointSize = u_sizeStart  + amount * (u_sizeEnd - u_sizeStart);
}
)";

#endif //FLUIDSOLVER_GRIDSCAL_VERT_H
