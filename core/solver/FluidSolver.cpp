//
// Created by austin on 2/28/16.
//

#include "FluidSolver.h"

float Fluid::Solver::g = -9.80665f;

Fluid::Solver::Solver(float particleSep) : particle_radius(particleSep) {

}

Fluid::Solver::~Solver() {
    delete _container;
    for (Particle* p : _particles) {
        delete p;
    }
}

void Fluid::Solver::setContainer(GeoObject *container) {
    _container = container;
}

void Fluid::Solver::addFluid(GeoObject *fluid) {
    Bound& b = fluid->bound();
    for (float x = b.minX(); x < b.maxX(); x += particle_radius) {
        for (float y = b.minY(); y < b.maxY(); y += particle_radius) {
            for (float z = b.minZ(); z < b.maxZ(); z += particle_radius) {
                glm::vec3 pos = glm::vec3(x, y, z);
                if (fluid->contains(pos)) {
                    Particle* p = new Particle();
                    p->pos = pos;
                    _particles.push_back(p);
                }
            }
        }
    }
}

void Fluid::Solver::update(float step) {
    presolve(step);
    solve(step);
    postsolve(step);
}

void Fluid::Solver::presolve(float step) {

}

void Fluid::Solver::solve(float step) {
    for (Particle* p : _particles) {
        p->vel += g * step;
        p->pos += p->vel * step;
    }
}

void Fluid::Solver::postsolve(float step) {

}
