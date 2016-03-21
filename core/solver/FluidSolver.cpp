//
// Created by austin on 2/28/16.
//

#include <iostream>
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

void FluidSolver::init() {
    glm::vec3 size = _container->bound().dim();
    glm::vec3 origin = (_container->bound().center() - size) / 2.f;
    _cell_size = 4*particle_radius;
    _MAC = MACGrid<std::vector<FluidParticle*> >(origin, size, _cell_size);

    for (FluidParticle &particle : _particles) {
        particle.cell = _MAC.indexOf(particle.pos);
        _MAC(particle.cell).push_back(&particle);
    }

    std::size_t velOffset = offsetof(FluidParticle, vel);
    std::size_t U_offset = velOffset + offsetof(glm::vec3, x);
    std::size_t V_offset = velOffset + offsetof(glm::vec3, y);
    std::size_t W_offset = velOffset + offsetof(glm::vec3, z);
    particleAttributeToGrid(U_offset, _MAC._gU_old, _cell_size, 0.f);
    particleAttributeToGrid(V_offset, _MAC._gV_old, _cell_size, 0.f);
    particleAttributeToGrid(W_offset, _MAC._gW_old, _cell_size, 0.f);
    particleAttributeToGrid(U_offset, _MAC._gU, _cell_size, 0.f);
    particleAttributeToGrid(V_offset, _MAC._gV, _cell_size, 0.f);
    particleAttributeToGrid(W_offset, _MAC._gW, _cell_size, 0.f);
}

void FluidSolver::projectVelocitiesToGrid() {
    std::size_t velOffset = offsetof(FluidParticle, vel);
    std::size_t U_offset = velOffset + offsetof(glm::vec3, x);
    std::size_t V_offset = velOffset + offsetof(glm::vec3, y);
    std::size_t W_offset = velOffset + offsetof(glm::vec3, z);

    _MAC._gU_old = _MAC._gU;
    _MAC._gV_old = _MAC._gV;
    _MAC._gW_old = _MAC._gW;
    particleAttributeToGrid(U_offset, _MAC._gU, _cell_size, 0.f);
    particleAttributeToGrid(V_offset, _MAC._gV, _cell_size, 0.f);
    particleAttributeToGrid(W_offset, _MAC._gW, _cell_size, 0.f);
}

void FluidSolver::transferVelocitiesToParticles() {
    float smooth = 1.f;
    for (FluidParticle &particle : _particles) {
        float vel = interpolateAttribute(particle.pos, _MAC._gU);
        float oldVel = interpolateAttribute(particle.pos, _MAC._gU_old);
        particle.vel.x = vel*smooth + (particle.vel.x +(vel - oldVel))*(1.f-smooth);
    }
    for (FluidParticle &particle : _particles) {
        float vel = interpolateAttribute(particle.pos, _MAC._gV);
        float oldVel = interpolateAttribute(particle.pos, _MAC._gV_old);
        particle.vel.y = vel*smooth + (particle.vel.y +(vel - oldVel))*(1.f-smooth);
    }
    for (FluidParticle &particle : _particles) {
        float vel = interpolateAttribute(particle.pos, _MAC._gW);
        float oldVel = interpolateAttribute(particle.pos, _MAC._gW_old);
        particle.vel.z = vel*smooth + (particle.vel.z +(vel - oldVel))*(1.f-smooth);
    }
}

void FluidSolver::gravitySolve(float step) {
    _MAC._gV.iterate([&](size_t i, size_t j, size_t k) {
        _MAC._gV(i, j, k) += g * step;
    });
}

void FluidSolver::updateParticlePositions(float step) {
    for (FluidParticle &particle : _particles) {
        particle.pos_old = particle.pos;
        particle.pos += particle.vel * step;
    }
}

void FluidSolver::resolveCollisions() {
    for (FluidParticle &particle : _particles) {
        glm::vec3 normal;
        if (_container->collides(particle.pos_old, particle.pos, normal)) {
            particle.col = glm::vec3(1,0,0);
            glm::vec3 mask = glm::vec3(1,1,1) - normal;
            particle.vel *= mask;
            particle.pos = particle.pos_old;
        }
    }
}

void FluidSolver::updateCells() {
    _MAC.iterate([&](size_t i, size_t j, size_t k) {
        _MAC(i,j,k).clear();
    });
    for (FluidParticle &particle : _particles) {
        particle.cell = _MAC.indexOf(particle.pos);
        if (particle.cell.x < 0 || particle.cell.y < 0 || particle.cell.z < 0) {
            std::cerr << "particle out of bounds" << std::endl;
        } else {
            _MAC(particle.cell).push_back(&particle);
        }
    }
}

template<class T> void FluidSolver::particleAttributeToGrid(std::size_t offset, Grid<T> &grid, float radius, T zeroVal) {
    std::size_t attributeSize = sizeof(T);
    std::size_t cellRadius = (size_t) glm::ceil(radius / _cell_size);

    grid.iterate([&](size_t I, size_t J, size_t K) {
        glm::vec3 gridPos = grid.positionOf(I,J,K);

        float totalDist = 0;
        for (size_t i = I - cellRadius; i < I + cellRadius; i++) {
            for (size_t j = J - cellRadius; j < J + cellRadius; j++) {
                for (size_t k = K - cellRadius; k < K + cellRadius; k++) {
                    for (FluidParticle* particle : _MAC(i,j,k)) {
                        float dist = glm::distance(particle->pos, gridPos);
                        if (dist < radius) {
                            totalDist += dist;
                        }
                    }
                }
            }
        }

        T temp;
        T gridVal = zeroVal;
        for (size_t i = I - cellRadius; i < I + cellRadius; i++) {
            for (size_t j = J - cellRadius; j < J + cellRadius; j++) {
                for (size_t k = K - cellRadius; k < K + cellRadius; k++) {
                    for (FluidParticle* particle : _MAC(i,j,k)) {
                        float dist = glm::distance(particle->pos, gridPos);
                        if (dist < radius) {
                            void* address = (void*)particle + offset;
                            std::memcpy(&temp, address, attributeSize);
                            gridVal += temp * dist / totalDist;
                        }
                    }
                }
            }
        }

        grid(I,J,K) = gridVal;
    });
}

template<class T> T FluidSolver::interpolateAttribute(const glm::vec3 &pos, Grid<T> &grid) {
    glm::vec3 idx = grid.fractionalIndexOf(pos);
    size_t i = (size_t) floor(idx.x);
    size_t j = (size_t) floor(idx.y);
    size_t k = (size_t) floor(idx.z);
    size_t I = (size_t) ceil(idx.x);
    size_t J = (size_t) ceil(idx.y);
    size_t K = (size_t) ceil(idx.z);
    
    T k1 = (idx.z-k) / (K-k) * grid(i,j,k) + (K-idx.z) / (K-k) * grid(i,j,K);
    T k2 = (idx.z-k) / (K-k) * grid(i,J,k) + (K-idx.z) / (K-k) * grid(i,J,K);
    T k3 = (idx.z-k) / (K-k) * grid(I,j,k) + (K-idx.z) / (K-k) * grid(I,j,K);
    T k4 = (idx.z-k) / (K-k) * grid(I,J,k) + (K-idx.z) / (K-k) * grid(I,J,K);

    T j1 = (idx.y-j) / (J-j) * k1 + (J-idx.y) / (J-j) * k2;
    T j2 = (idx.y-j) / (J-j) * k3 + (J-idx.y) / (J-j) * k4;

    return (idx.x-i) / (I-i) * j1 + (I-idx.x) / (I-i) * j2;
}

void FluidSolver::update(float step) {
    projectVelocitiesToGrid();
    gravitySolve(step);
    transferVelocitiesToParticles();
    updateParticlePositions(step);
    resolveCollisions();
    updateCells();
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