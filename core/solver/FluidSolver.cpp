//
// Created by austin on 2/28/16.
//

#include "FluidSolver.h"

float FluidSolver::g = -9.80665f;

FluidSolver::FluidSolver(float particleSep) : particle_radius(particleSep) {
    FluidParticle p1;
    p1.pos = glm::vec3(-0.2, -0.2, 0);

    FluidParticle p2;
    p2.pos = glm::vec3(0, 0.8, 0);

    FluidParticle p3;
    p3.pos = glm::vec3(0.8, 0.8, 0);

    _particles.push_back(p1);
    _particles.push_back(p2);
    _particles.push_back(p3);
}

FluidSolver::~FluidSolver() {
    delete _container;
}

void FluidSolver::setContainer(GeoObject *container) {
    _container = container;
}

void FluidSolver::addFluid(GeoObject *fluid) {
    Bound& b = fluid->bound();
    for (float x = b.minX(); x < b.maxX(); x += particle_radius) {
        for (float y = b.minY(); y < b.maxY(); y += particle_radius) {
            for (float z = b.minZ(); z < b.maxZ(); z += particle_radius) {
                glm::vec3 pos = glm::vec3(x, y, z);
                if (fluid->contains(pos)) {
                    FluidParticle p;
                    p.pos = pos;
                    _particles.push_back(p);
                }
            }
        }
    }
}

void FluidSolver::update(float step) {
    presolve(step);
    solve(step);
    postsolve(step);
}

void FluidSolver::presolve(float step) {

}

void FluidSolver::solve(float step) {
    for (FluidParticle &p : _particles) {
        p.vel += g * step;
        p.pos += p.vel * step;
    }
}

void FluidSolver::postsolve(float step) {

}
