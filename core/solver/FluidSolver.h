//
// Created by austin on 2/28/16.
//

#ifndef FLUIDSOLVER_FLUIDSOLVER_H
#define FLUIDSOLVER_FLUIDSOLVER_H

#include <core/util/math.h>
#include <core/geometry/GeoObject.h>
#include <vector>

struct FluidParticle {
    glm::vec3 pos;
    glm::vec3 vel;
    glm::vec3 col;

    FluidParticle() {
        pos = glm::vec3(0);
        vel = glm::vec3(0);
        col = glm::vec3(0.5f, 0.5f, 1.f);
    }
};

class FluidSolver {
    friend class ParticlesPainter;
public:
    FluidSolver(float particleSep);
    ~FluidSolver();

    void setContainer(GeoObject* container);
    void addFluid(GeoObject* fluid);

    void update(float step = 0.04166f);

    GeoObject* _container;

private:
    std::vector<FluidParticle> _particles;
    float particle_radius;

    void presolve(float step);
    void solve(float step);
    void postsolve(float step);

    static float g;
};


#endif //FLUIDSOLVER_FLUIDSOLVER_H
