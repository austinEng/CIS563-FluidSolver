//
// Created by austin on 3/22/16.
//

#ifndef FLUIDSOLVER_GRIDATTR_GEO_H
#define FLUIDSOLVER_GRIDATTR_GEO_H

const char* gridAttr_geo = R"(
#version 150

layout(points) in;
layout(line_strip, max_vertices = 2) out;

uniform mat4 u_viewProj;
uniform vec3 u_vec;

in float g_data[];
out float f_scale;

void main() {
    f_scale = abs(g_data[0]);

    gl_Position = u_viewProj * gl_in[0].gl_Position;
    EmitVertex();

    gl_Position = u_viewProj * (gl_in[0].gl_Position + g_data[0]*vec4(u_vec, 0));
    EmitVertex();
    EndPrimitive();
}
)";

#endif //FLUIDSOLVER_GRIDATTR_GEO_H
