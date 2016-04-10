//
// Created by austin on 2/29/16.
//

#include "Camera.h"

Camera::Camera(int w, int h) :
        zoom(25),
        eye(glm::vec3(0,-5,zoom)),
        tgt(glm::vec3(0,-5,0)),
        width(w),
        height(h),
        fovy(45),
        near_clip(0.001f),
        far_clip(1000.f),
        world_up(glm::vec3(0,1,0)),
        look(tgt - eye),
        right(glm::cross(look, world_up)),
        up(glm::cross(right, look)),
        rotation(glm::mat4(1.f)) {

    rotation = glm::rotate(rotation, -PI/6, glm::vec3(1,0,0));
    rotation = glm::rotate(rotation, -PI/8, glm::vec3(0,1,0));

    recomputeEye();
    recompute();
}

glm::mat4 Camera::viewProj() {
    return _viewProj;
}

void Camera::recompute() {
    float aspect = (float)width/height;
    look = glm::normalize(tgt - eye);
    right = glm::cross(look, world_up);
    up = glm::cross(right, look);

    _viewProj = glm::perspective(fovy, aspect, near_clip, far_clip) * glm::lookAt(eye, tgt, up);
}

void Camera::resize(int w, int h) {
    width = w;
    height = h;
    recompute();
}

void Camera::recomputeEye() {
    eye = glm::vec3(rotation * glm::vec4(0,0,zoom,0)) + tgt;
    look = glm::normalize(tgt - eye);
    right = glm::normalize(glm::cross(look, up));
    up = glm::normalize(glm::cross(right, look));
}
