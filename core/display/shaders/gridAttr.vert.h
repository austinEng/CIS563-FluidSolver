//
// Created by austin on 2/28/16.
//

#ifndef FLUIDSOLVER_GRIDATTR_VERT_H
#define FLUIDSOLVER_GRIDATTR_VERT_H

const char* gridAttr_vert = R"(
#version 150

uniform float u_cellSize;
uniform ivec3 u_cellCount;
uniform vec3 u_origin;

in float v_data;
out float g_data;

void main() {
    g_data = v_data;

    int i = int(mod(gl_VertexID, u_cellCount.z));
    int j = int(mod(gl_VertexID / u_cellCount.z, u_cellCount.y));
    int k = int(gl_VertexID / (u_cellCount.y * u_cellCount.z));

    vec3 pos = vec3(float(i), float(j), float(k)) * u_cellSize + u_origin;

    gl_Position = vec4(pos, 1);
}
)";

#endif