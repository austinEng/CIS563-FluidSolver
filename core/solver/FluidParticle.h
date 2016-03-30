//
// Created by austin on 3/20/16.
//

#ifndef FLUIDSOLVER_FLUIDPARTICLE_H
#define FLUIDSOLVER_FLUIDPARTICLE_H

#include <core/util/math.h>

struct FluidParticle {
    glm::vec3 pos;
    glm::vec3 pos_old;
    glm::vec3 vel;
    glm::vec3 col;
    glm::ivec3 cell;

    FluidParticle() {
        pos = glm::vec3(0);
        pos_old = glm::vec3(0);
        vel = glm::vec3(0,0,0);
        col = glm::vec3(0.5f, 0.5f, 1.f);
    }
};

#endif //FLUIDSOLVER_FLUIDPARTICLE_H
