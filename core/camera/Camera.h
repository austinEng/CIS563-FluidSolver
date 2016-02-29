//
// Created by austin on 2/29/16.
//

#ifndef FLUIDSOLVER_CAMERA_H
#define FLUIDSOLVER_CAMERA_H

#include <core/util/math.h>

class Camera {
private:
    glm::vec3 world_up;
public:
    float zoom;
    glm::vec3 eye;
    glm::vec3 tgt;

    Camera(int w, int h);
    glm::mat4 viewProj();
    void recomputeEye();
    void recompute();
    void resize(int w, int h);

    glm::vec3 look, up, right;
    int width, height;
    float fovy, near_clip, far_clip;
    glm::mat4 rotation;

private:

    glm::mat4 _viewProj;
};


#endif //FLUIDSOLVER_CAMERA_H
