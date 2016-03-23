//
// Created by austin on 2/28/16.
//

#include <iostream>
#include "FluidSolver.h"
#include <core/util/hacks.h>
#include <tbb/parallel_for.h>
#include <core/util/flags.h>
#include <tbb/parallel_invoke.h>
#include <tbb/parallel_reduce.h>

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
                glm::vec3 pos = glm::vec3(x, y, z) + glm::vec3(_cell_size)/2.f;
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

#ifdef USETBB
    tbb::parallel_invoke(
            [&]() {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, _particles.size()),
                                  [&](const tbb::blocked_range<size_t> &r) {
                                      for (size_t i = r.begin(); i != r.end(); ++i) {
                                          FluidParticle &particle = _particles[i];
                                          float vel = interpolateAttribute(particle.pos, _MAC._gU);
                                          float oldVel = interpolateAttribute(particle.pos, _MAC._gU_old);
                                          particle.vel.x = vel*smooth + (particle.vel.x +(vel - oldVel))*(1.f-smooth);
                                      }
                                  });
            },
            [&]() {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, _particles.size()),
                                  [&](const tbb::blocked_range<size_t> &r) {
                                      for (size_t i = r.begin(); i != r.end(); ++i) {
                                          FluidParticle &particle = _particles[i];
                                          float vel = interpolateAttribute(particle.pos, _MAC._gV);
                                          float oldVel = interpolateAttribute(particle.pos, _MAC._gV_old);
                                          particle.vel.y = vel*smooth + (particle.vel.y +(vel - oldVel))*(1.f-smooth);
                                      }
                                  });
            },
            [&]() {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, _particles.size()),
                                  [&](const tbb::blocked_range<size_t> &r) {
                                      for (size_t i = r.begin(); i != r.end(); ++i) {
                                          FluidParticle &particle = _particles[i];
                                          float vel = interpolateAttribute(particle.pos, _MAC._gW);
                                          float oldVel = interpolateAttribute(particle.pos, _MAC._gW_old);
                                          particle.vel.z = vel*smooth + (particle.vel.z +(vel - oldVel))*(1.f-smooth);
                                      }
                                  });
            }
    );
#else
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
#endif
}

void FluidSolver::gravitySolve(float step) {

#ifdef USETBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, _particles.size()), [&](const tbb::blocked_range<size_t> &r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            FluidParticle &particle = _particles[i];
            particle.vel.y += g*step;
        }
    });
#else
    for (FluidParticle &particle : _particles) {
        particle.vel.y += g*step;
    }
#endif

}

void FluidSolver::updateParticlePositions(float step) {

#ifdef USETBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, _particles.size()), [&](const tbb::blocked_range<size_t> &r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            FluidParticle &particle = _particles[i];
            particle.pos_old = particle.pos;
            particle.pos += particle.vel * step;
        }
    });
#else
    for (FluidParticle &particle : _particles) {
        particle.pos_old = particle.pos;
        particle.pos += particle.vel * step;
    }
#endif
}

void FluidSolver::resolveCollisions() {

#ifdef USETBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, _particles.size()), [&](const tbb::blocked_range<size_t> &r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            FluidParticle &particle = _particles[i];
            glm::vec3 normal;
            if (_container->collides(particle.pos_old, particle.pos, normal)) {
                particle.col = glm::vec3(1,0,0);
                glm::vec3 mask = glm::vec3(1,1,1) - glm::abs(normal);
                particle.vel *= mask;
                particle.pos = particle.pos_old;
            }
        }
    });
#else
    for (FluidParticle &particle : _particles) {
        glm::vec3 normal;
        if (_container->collides(particle.pos_old, particle.pos, normal)) {
            particle.col = glm::vec3(1,0,0);
            glm::vec3 mask = glm::vec3(1,1,1) - glm::abs(normal);
            particle.vel *= mask;
            particle.pos = particle.pos_old;
        }
    }
#endif
}

void FluidSolver::updateCells() {

    for (FluidParticle &particle : _particles) {
        particle.cell = _MAC.indexOf(particle.pos);
        if (_MAC.checkIdx(particle.cell)) {
            _MAC(particle.cell).push_back(&particle);
        } else {
            //std::cerr << "particle out of bounds" << std::endl;
        }

    }

}

template<class T> void FluidSolver::particleAttributeToGrid(std::size_t offset, Grid<T> &grid, float radius, T zeroVal) {
    std::size_t attributeSize = sizeof(T);
    std::size_t cellRadius = (size_t) glm::ceil(radius / _cell_size);

    /*grid.clear(zeroVal);
    Grid<float> distGrid(grid);

    iterParticles([&](FluidParticle &particle) {
        size_t I,J,K;
        grid.indexOf(particle.pos, I, J, K);
        glm::vec3 gridPos = distGrid.positionOf(I,J,K);

        distGrid.iterateNeighborhood(I,J,K,cellRadius, [&](size_t i, size_t j, size_t k) {
            float dist = glm::distance(particle.pos, gridPos);
            distGrid(i,j,k) += dist * (dist < radius);
        }, false);
    }, false);

    iterParticles([&](FluidParticle &particle) {
        size_t I,J,K;
        grid.indexOf(particle.pos, I, J, K);
        glm::vec3 gridPos = distGrid.positionOf(I,J,K);

        distGrid.iterateNeighborhood(I,J,K,cellRadius, [&](size_t i, size_t j, size_t k) {
            float dist = glm::distance(particle.pos, gridPos);
            T temp;
            void *address = (void *) &particle + offset;
            std::memcpy(&temp, address, attributeSize);
            grid(i,j,k) += temp * (dist / distGrid(i,j,k)) * (dist < radius);
        }, false);
    }, false);*/

    grid.iterate([&](size_t I, size_t J, size_t K) {

        glm::vec3 gridPos = grid.positionOf(I,J,K);

        size_t si, ei, sj, ej, sk, ek;
        _MAC.getNeighboorhood(I, J, K, cellRadius, si, ei, sj, ej, sk, ek);

        float totalDist = 0.f;
        for (size_t i = si; i < ei; i++) {
            for (size_t j = sj; j < ej; j++) {
                for (size_t k = sk; k < ek; k++) {
#ifdef USETBB
                    typedef std::vector<FluidParticle*>::iterator range_iterator;
                    totalDist += tbb::parallel_reduce(
                            tbb::blocked_range<range_iterator>(_MAC(i,j,k).begin(), _MAC(i,j,k).end(), 20), 0.f,
                            [&](const tbb::blocked_range<range_iterator> &r, float init)->T {
                                for (range_iterator p = r.begin(); p != r.end(); p++) {
                                    float dist = glm::distance((*p)->pos, gridPos);
                                    init += dist * (dist < radius);
                                }
                                return init;
                            }, std::plus<float>());
#else
                    for (FluidParticle* particle : _MAC(i,j,k)) {
                        float dist = glm::distance(particle->pos, gridPos);
                        if (dist < radius) {
                            totalDist += dist;
                        }
                    }
#endif
                }
            }
        }

        if (totalDist == 0) {
            grid(I,J,K) = zeroVal;
            return;
        }

        T temp;
        T gridVal = zeroVal;
        for (size_t i = si; i < ei; i++) {
            for (size_t j = sj; j < ej; j++) {
                for (size_t k = sk; k < ek; k++) {
                    for (FluidParticle *particle : _MAC(i,j,k)) {
                        float dist = glm::distance(particle->pos, gridPos);
                        void *address = (void *) particle + offset;
                        std::memcpy(&temp, address, attributeSize);
                        gridVal += temp * (dist / totalDist) * (dist < radius);
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

    // this is reverse from what is expected because we want smaller value (closer distance) to have larger influence
    k1 = MATHIFELSE((K-idx.z) * grid(i,j,k) + (idx.z-k) * grid(i,j,K), grid(i,j,k), k==K);
    k2 = MATHIFELSE((K-idx.z) * grid(i,J,k) + (idx.z-k) * grid(i,J,K), grid(i,J,k), k==K);
    k3 = MATHIFELSE((K-idx.z) * grid(I,j,k) + (idx.z-k) * grid(I,j,K), grid(I,j,k), k==K);
    k4 = MATHIFELSE((K-idx.z) * grid(I,J,k) + (idx.z-k) * grid(i,J,K), grid(I,J,k), k==K);

    j1 = MATHIFELSE((J-idx.y) * k1 + (idx.y-j) * k2, k1, j==J);
    j2 = MATHIFELSE((J-idx.y) * k3 + (idx.y-j) * k4, k3, j==J);

    val = MATHIFELSE((I-idx.x) * j1 + (idx.x-i) * j2, j1, i==I);

    return val;
}

void FluidSolver::update(float step) {
    projectVelocitiesToGrid();
    // pressure solve
    transferVelocitiesToParticles();
    gravitySolve(step);
    updateParticlePositions(step);
    resolveCollisions();
    updateCells();

    frame++;
}

void FluidSolver::iterParticles(const std::function<void(FluidParticle &particle)> &cb, bool parallel) {
#ifdef USETBB
    if (parallel) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, _particles.size()), [&](const tbb::blocked_range<size_t> &r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                FluidParticle &particle = _particles[i];
                cb(particle);
            }
        });
    } else {
        for (FluidParticle &particle : _particles) {
            cb(particle);
        }
    }
#else
    for (FluidParticle &particle : _particles) {
        cb(particle);
    }
#endif
}