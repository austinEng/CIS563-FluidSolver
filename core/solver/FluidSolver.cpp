//
// Created by austin on 2/28/16.
//

#include <iostream>
#include "FluidSolver.h"
#include <core/util/hacks.h>

float FluidSolver::g = -9.80665f;

FluidSolver::FluidSolver(float particleSep, float gridSize) : particle_radius(particleSep), _cell_size(gridSize), frame(0) {
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
        //std::cout << particle.vel.y << std::endl;
    }
    for (FluidParticle &particle : _particles) {
        float vel = interpolateAttribute(particle.pos, _MAC._gW);
        float oldVel = interpolateAttribute(particle.pos, _MAC._gW_old);
        particle.vel.z = vel*smooth + (particle.vel.z +(vel - oldVel))*(1.f-smooth);
    }
}

void FluidSolver::gravitySolve(float step) {
    for (FluidParticle &particle : _particles) {
        particle.vel.y += g*step;
    }
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
            //std::cout << glm::to_string(particle.pos_old) << "; " << glm::to_string(particle.pos) << std::endl;
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
        size_t startI = MATHIFELSE(I - cellRadius, I, I == 0);
        size_t startJ = MATHIFELSE(J - cellRadius, J, J == 0);
        size_t startK = MATHIFELSE(K - cellRadius, K, K == 0);
        for (size_t i = startI; i <= I + cellRadius; i++) {
            for (size_t j = startJ; j <= J + cellRadius; j++) {
                for (size_t k = startK; k <= K + cellRadius; k++) {
                    if (_MAC.checkIdx(i,j,k)) {
                        for (FluidParticle* particle : _MAC(i,j,k)) {
                            float dist = glm::distance(particle->pos, gridPos);
                            if (dist < radius) {
                                totalDist += dist;
                            }
                        }
                    }
                }
            }
        }


        T temp;
        T gridVal = zeroVal;
        for (size_t i = I - cellRadius; i <= I + cellRadius; i++) {
            for (size_t j = J - cellRadius; j <= J + cellRadius; j++) {
                for (size_t k = K - cellRadius; k <= K + cellRadius; k++) {
                    if (_MAC.checkIdx(i,j,k)) {
                        for (FluidParticle *particle : _MAC(i, j, k)) {
                            float dist = glm::distance(particle->pos, gridPos);
                            if (dist < radius) {
                                void *address = (void *) particle + offset;
                                std::memcpy(&temp, address, attributeSize);
                                gridVal += temp * dist / totalDist;
                            }
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

    /*float divisor = 1.f;
    if (K-k != 0.f) divisor = K-k;
    T k1 = (idx.z-k) / divisor * grid(i,j,k) + (K-idx.z) / divisor * grid(i,j,K);
    T k2 = (idx.z-k) / divisor * grid(i,J,k) + (K-idx.z) / divisor * grid(i,J,K);
    T k3 = (idx.z-k) / divisor * grid(I,j,k) + (K-idx.z) / divisor * grid(I,j,K);
    T k4 = (idx.z-k) / divisor * grid(I,J,k) + (K-idx.z) / divisor * grid(I,J,K);

    divisor = 1.f;
    if (J-j != 0.f) divisor = J-j;
    T j1 = (idx.y-j) / divisor * k1 + (J-idx.y) / divisor * k2;
    T j2 = (idx.y-j) / divisor * k3 + (J-idx.y) / divisor * k4;

    divisor = 1.f;
    if (I-i != 0.f) divisor = I-i;
    T val = (idx.x-i) / divisor * j1 + (I-idx.x) / divisor * j2;*/

    T k1, k2, k3, k4, j1, j2, val;

    /*if (k == K) {
        k1 = grid(i,j,k);
        k2 = grid(i,J,k);
        k3 = grid(I,j,k);
        k4 = grid(I,J,k);
    } else {
        k1 = (idx.z-k) * grid(i,j,k) + (K-idx.z) * grid(i,j,K);
        k2 = (idx.z-k) * grid(i,J,k) + (K-idx.z) * grid(i,J,K);
        k3 = (idx.z-k) * grid(I,j,k) + (K-idx.z) * grid(I,j,K);
        k4 = (idx.z-k) * grid(I,J,k) + (K-idx.z) * grid(I,J,K);
    }

    if (j == J) {
        j1 = k1;
        j2 = k3;
    } else {
        j1 = (idx.y-j) * k1 + (J-idx.y) * k2;
        j2 = (idx.y-j) * k3 + (J-idx.y) * k4;
    }

    if (i == I) {
        val = j1;
    } else {
        val = (idx.x-i) * j1 + (I-idx.x) * j2;
    }*/

    k1 = MATHIFELSE((idx.z-k) * grid(i,j,k) + (K-idx.z) * grid(i,j,K), grid(i,j,k), k==K);
    k2 = MATHIFELSE((idx.z-k) * grid(i,J,k) + (K-idx.z) * grid(i,J,K), grid(i,J,k), k==K);
    k3 = MATHIFELSE((idx.z-k) * grid(I,j,k) + (K-idx.z) * grid(I,j,K), grid(I,j,k), k==K);
    k4 = MATHIFELSE((idx.z-k) * grid(I,J,k) + (K-idx.z) * grid(i,J,K), grid(I,J,k), k==K);

    j1 = MATHIFELSE((idx.y-j) * k1 + (J-idx.y) * k2, k1, j==J);
    j2 = MATHIFELSE((idx.y-j) * k3 + (J-idx.y) * k4, k3, j==J);

    val = MATHIFELSE((idx.x-i) * j1 + (I-idx.x) * j2, j1, i==I);

    return val;
}

void FluidSolver::update(float step) {
    /*std::cout << "================PARTICLE COUNT===============" << std::endl;
    _MAC.iterate([&](size_t i, size_t j, size_t k) {
        std::cout << i << "," << j << "," << k << "; " << _MAC(i,j,k).size() << std::endl;
    });

    projectVelocitiesToGrid();

    std::cout << "================U VELOCITY===============" << std::endl;
    _MAC._gU.iterate([&](size_t i, size_t j, size_t k) {
        std::cout << i << "," << j << "," << k << "; " << _MAC._gU(i,j,k) << std::endl;
    });
    std::cout << "================V VELOCITY===============" << std::endl;
    _MAC._gV.iterate([&](size_t i, size_t j, size_t k) {
        std::cout << i << "," << j << "," << k << "; " << _MAC._gV(i,j,k) << std::endl;
    });
    std::cout << "================W VELOCITY===============" << std::endl;
    _MAC._gW.iterate([&](size_t i, size_t j, size_t k) {
        std::cout << i << "," << j << "," << k << "; " << _MAC._gW(i,j,k) << std::endl;
    });*/

    transferVelocitiesToParticles();

    /*std::cout << "================PARTICLE VELOCITY===============" << std::endl;
    for (FluidParticle &particle : _particles) {
        std::cout << glm::to_string(particle.pos) << "; " << glm::to_string(particle.vel) << std::endl;
    }*/

    gravitySolve(step);
    updateParticlePositions(step);
    resolveCollisions();
    updateCells();



    frame++;
}