//
// Created by austin on 2/28/16.
//

#include "FluidSolver.h"

float FluidSolver::g = -9.80665f;

FluidSolver::FluidSolver(float particleSep) : particle_radius(particleSep) {
}

FluidSolver::~FluidSolver() {
    delete _container;
}

void FluidSolver::setContainer(GeoObject *container) {
    _container = container;
}

/*
 * Loop over fluid bounds to generate particles
 */
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
    std::vector<FluidParticle> tempParticles(_particles.size());
    for (unsigned int i = 0; i < _particles.size(); i++) {
        tempParticles[i].vel = _particles[i].vel + glm::vec3(0, g * step, 0);
        tempParticles[i].pos = _particles[i].pos + _particles[i].vel * step;

        if (_container->collides(_particles[i].pos, tempParticles[i].pos)) {
            _particles[i].col = glm::vec3(1,0,0);
        }
        tempParticles[i].col = _particles[i].col;
    }
    swap(_particles, tempParticles);
}

void FluidSolver::postsolve(float step) {

}
