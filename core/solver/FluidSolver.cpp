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
#include <tbb/blocked_range3d.h>
#include <tbb/partitioner.h>

float FluidSolver::g = -9.80665f;

FluidSolver::FluidSolver(float particleSep, float gridSize) : particle_radius(particleSep), _cell_size(gridSize), frame(0) {
}

FluidSolver::~FluidSolver() {
    delete _container;
}

void FluidSolver::setContainer(GeoObject *container) {
    _container = container;
    glm::vec3 size = _container->bound().dim();
    glm::vec3 origin = (_container->bound().center() - size) / 2.f;
    _MAC = MACGrid<std::vector<FluidParticle*> >(
            origin - glm::vec3(_cell_size, _cell_size, _cell_size),
            size + 2.f*glm::vec3(_cell_size, _cell_size, _cell_size),
            _cell_size
    );

    std::function<void(size_t, size_t, size_t)> setSolid = [&](size_t i, size_t j, size_t k) {
        _MAC._gType(i,j,k) = SOLID;
    };

    _MAC._gType.iterateRegion(0,0,0, 1,_MAC._gType.countY(),_MAC._gType.countZ(), setSolid);
    _MAC._gType.iterateRegion(_MAC._gType.countX()-1,0,0, _MAC._gType.countX(),_MAC._gType.countY(),_MAC._gType.countZ(), setSolid);
    _MAC._gType.iterateRegion(0,0,0, _MAC._gType.countX(),1,_MAC._gType.countZ(), setSolid);
    _MAC._gType.iterateRegion(0,_MAC._gType.countY()-1,0, _MAC._gType.countX(),_MAC._gType.countY(),_MAC._gType.countZ(), setSolid);
    _MAC._gType.iterateRegion(0,0,0, _MAC._gType.countX(),_MAC._gType.countY(),1, setSolid);
    _MAC._gType.iterateRegion(0,0,_MAC._gType.countZ()-1, _MAC._gType.countX(),_MAC._gType.countY(),_MAC._gType.countZ(), setSolid);

    /*_MAC = MACGrid<std::vector<FluidParticle*> >(
            origin,
            size,
            _cell_size
    );*/
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
                    _MAC._gType.at(pos) = FLUID;
                }
            }
        }
    }
}

void FluidSolver::init() {

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
    float smooth = 0.05f;

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

void FluidSolver::enforceBoundary() {
    _MAC._gType.iterate([&](size_t i, size_t j, size_t k) {
        switch (_MAC._gType(i,j,k)) {
            case EMPTY:break;
            case FLUID:break;
            case SOLID:
                size_t si, ei, sj, ej, sk, ek;
                _MAC._gType.getNeighboorhood(i,j,k,1,si,ei,sj,ej,sk,ek);
                _MAC._gU(i,j,k) = std::min(_MAC._gU(i,j,k), 0.f);
                _MAC._gU(ei-1,j,k) = std::max(_MAC._gU(ei-1,j,k), 0.f);
                _MAC._gV(i,j,k) = std::min(_MAC._gV(i,j,k), 0.f);
                _MAC._gV(i,ej-1,k) = std::max(_MAC._gV(i,ej-1,k), 0.f);
                _MAC._gW(i,j,k) = std::min(_MAC._gW(i,j,k), 0.f);
                _MAC._gW(i,j,ek-1) = std::max(_MAC._gW(i,j,ek-1), 0.f);
                break;
            default:break;
        }
    });
}

void FluidSolver::gravitySolve(float step) {
    _MAC._gV.iterate([&](size_t i, size_t j, size_t k) {
        _MAC._gV(i,j,k) += g*step;
    });
}

void FluidSolver::extrapolateVelocity() {
    _MAC._gType.iterate([&](size_t i, size_t j, size_t k) {
        if (_MAC._gType(i,j,k) != FLUID) {
            if (i > 0 && _MAC._gType(i-1,j,k) == FLUID) {

            }
            if (i+1 < _MAC._gType.countX() && _MAC._gType(i+1,j,k) == FLUID) {

            }
            if (j > 0 && _MAC._gType(i,j-1,k) == FLUID) {

            }
            if (j+1 < _MAC._gType.countY() && _MAC._gType(i,j+1,k) == FLUID) {

            }
            if (k > 0 && _MAC._gType(i,j,k-1) == FLUID) {

            }
            if (k+1 < _MAC._gType.countZ() && _MAC._gType(i,j,k+1) == FLUID) {

            }

            //size_t si, ei, sj, ej, sk, ek;
            int count = 0;
            //_MAC._gType.getNeighboorhood(i,j,k,1,si,ei,sj,ej,sk,ek) {


            //_MAC._gU(i,j,k) = 2.f;
        }
    });
}

void FluidSolver::updateParticlePositions(float step) {

#ifdef USETBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, _particles.size()), [&](const tbb::blocked_range<size_t> &r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            FluidParticle &particle = _particles[i];
            particle.pos_old = particle.pos;

            glm::vec3 k1 = step*particle.vel; // step * velAt(pos)
            particle.pos += step * glm::vec3(
                    interpolateAttribute(particle.pos + 0.5f*k1, _MAC._gU),
                    interpolateAttribute(particle.pos + 0.5f*k1, _MAC._gV),
                    interpolateAttribute(particle.pos + 0.5f*k1, _MAC._gW)
            );
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
                //particle.vel *= mask;
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
            //particle.vel *= mask;
            particle.pos = particle.pos_old;
        }
    }
#endif
}

void FluidSolver::updateCells() {
    _MAC.clear(std::vector<FluidParticle*>());
    _MAC._gType.clear(EMPTY);
    for (FluidParticle &particle : _particles) {
        particle.cell = _MAC.indexOf(particle.pos);
        _MAC._gType(particle.cell) = FLUID;
        if (_MAC.checkIdx(particle.cell)) {
            _MAC(particle.cell).push_back(&particle);
        } else {
            //std::cerr << "particle out of bounds" << std::endl;
        }

    }

    std::function<void(size_t, size_t, size_t)> setSolid = [&](size_t i, size_t j, size_t k) {
        _MAC._gType(i,j,k) = SOLID;
    };

    _MAC._gType.iterateRegion(0,0,0, 1,_MAC._gType.countY(),_MAC._gType.countZ(), setSolid);
    _MAC._gType.iterateRegion(_MAC._gType.countX()-1,0,0, _MAC._gType.countX(),_MAC._gType.countY(),_MAC._gType.countZ(), setSolid);
    _MAC._gType.iterateRegion(0,0,0, _MAC._gType.countX(),1,_MAC._gType.countZ(), setSolid);
    _MAC._gType.iterateRegion(0,_MAC._gType.countY()-1,0, _MAC._gType.countX(),_MAC._gType.countY(),_MAC._gType.countZ(), setSolid);
    _MAC._gType.iterateRegion(0,0,0, _MAC._gType.countX(),_MAC._gType.countY(),1, setSolid);
    _MAC._gType.iterateRegion(0,0,_MAC._gType.countZ()-1, _MAC._gType.countX(),_MAC._gType.countY(),_MAC._gType.countZ(), setSolid);

}

float kernel(float r, float h) {
    float e = r/h;
    return 1.f/(PI*h*h*h) * MATHIFELSE(
            1.f - 3.f/2.f * e*e  + 3.f/4.f * e*e*e,
            MATHIFELSE(
                    1.f/4.f * (2-e)*(2-e)*(2-e),
                    0,
                    e > 2
            ),
            e > 1
    );
}


template<typename T> void FluidSolver::particleAttributeToGrid(std::size_t offset, Grid<T> &grid, float radius, T zeroVal) {
    std::size_t attributeSize = sizeof(T);
    std::size_t cellRadius = (size_t) glm::ceil(radius / _cell_size);

#ifdef SPLATTING // if splatting
    grid.clear(zeroVal);
    std::vector<float> distances(grid.countX() * grid.countY() * grid.countZ());

    iterParticles([&](FluidParticle &particle) {
        size_t I,J,K;
        grid.indexOf(particle.pos, I, J, K);
        glm::vec3 gridPos = grid.positionOf(I,J,K);

        grid.iterateNeighborhood(I,J,K,cellRadius, [&](size_t i, size_t j, size_t k) {
            float dist = glm::distance(particle.pos, gridPos);
            size_t idx = grid.fromIJK(i,j,k);
            distances[idx] += dist * (dist < radius);
        });
    }, false);

    iterParticles([&](FluidParticle &particle) {
        size_t I,J,K;
        grid.indexOf(particle.pos, I, J, K);
        glm::vec3 gridPos = grid.positionOf(I,J,K);

        grid.iterateNeighborhood(I,J,K,cellRadius, [&](size_t i, size_t j, size_t k) {
            float dist = glm::distance(particle.pos, gridPos);
            size_t idx = grid.fromIJK(i,j,k);
            T temp;
            void *address = (void *) &particle + offset;
            std::memcpy(&temp, address, attributeSize);
            grid(i,j,k) += temp * (dist / distances[idx]) * (dist < radius);
        });
    }, false);

#else // else splatting

/*#ifdef USETBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, grid.size()), [&](const tbb::blocked_range<size_t> &r) {
        for (size_t idx = r.begin(); idx != r.end(); idx++) {
#else
    for (size_t idx = 0; idx < grid.size(); idx++) {
#endif
        size_t I, J, K;
        grid.toIJK(idx, I, J, K);
        glm::vec3 gridPos = grid.positionOf(I,J,K);

        size_t si, ei, sj, ej, sk, ek;
        _MAC.getNeighboorhood(I, J, K, cellRadius, si, ei, sj, ej, sk, ek);

        float totalWeight = 0.f;
        for (size_t i = si; i < ei; i++) {
            for (size_t j = sj; j < ej; j++) {
                for (size_t k = sk; k < ek; k++) {
                    for (FluidParticle *particle : _MAC(i,j,k)) {
                        float dist = glm::distance(particle->pos, gridPos);
                        float weight = kernel(dist, radius/2.f);
                        totalWeight += weight;
                    }
                }
            }
        }

        if (totalWeight == 0) {
            grid(I,J,K) = zeroVal;
            continue;
        }

        T temp;
        T gridVal = zeroVal;
        for (size_t i = si; i < ei; i++) {
            for (size_t j = sj; j < ej; j++) {
                for (size_t k = sk; k < ek; k++) {
                    for (FluidParticle *particle : _MAC(i,j,k)) {
                        float dist = glm::distance(particle->pos, gridPos);
                        float weight = kernel(dist, radius/2.f);
                        void *address = (void *) particle + offset;
                        std::memcpy(&temp, address, attributeSize);
                        gridVal += temp * (weight / totalWeight);
                    }
                }
            }
        }

        grid(I,J,K) = gridVal;
#ifdef USETBB
    });
#else
    }
#endif*/

    grid.iterate([&](size_t I, size_t J, size_t K) {
        glm::vec3 gridPos = grid.positionOf(I,J,K);

        size_t mI, mJ, mK, si, ei, sj, ej, sk, ek;
        _MAC.indexOf(gridPos, mI, mJ, mK);
        _MAC.getNeighboorhood(mI, mJ, mK, cellRadius, si, ei, sj, ej, sk, ek);

        float totalWeight = 0.f;
//#define TOOMUCH
#ifdef TOOMUCH
        totalWeight += tbb::parallel_reduce(
            tbb::blocked_range3d<size_t>(si,ei,sj,ej,sk,ek), 0.f,
            [&](const tbb::blocked_range3d<size_t> &r, float init)->float {
                for(size_t i=r.pages().begin(), i_end=r.pages().end(); i<i_end; i++){
                    for (size_t j=r.rows().begin(), j_end=r.rows().end(); j<j_end; j++){
                        for (size_t k=r.cols().begin(), k_end=r.cols().end(); k<k_end; k++){
                            typedef std::vector<FluidParticle*>::iterator range_iterator;
                            init += tbb::parallel_reduce(
                                    tbb::blocked_range<range_iterator>(_MAC(i,j,k).begin(), _MAC(i,j,k).end()), 0.f,
                                    [&](const tbb::blocked_range<range_iterator> &r2, float init2)->float {
                                        for (range_iterator p = r2.begin(); p != r2.end(); p++) {
                                            FluidParticle* particle = *p;
                                            float dist = glm::distance(particle->pos, gridPos);
                                            float weight = kernel(dist, radius/2.f);
                                            init2 += weight;
                                        }
                                        return init2;
                                    }, std::plus<float>()
                            );
                        }
                    }
                }
                return init;
            },
            std::plus<float>()
        );

        T gridVal = zeroVal;
        gridVal += tbb::parallel_reduce(
                tbb::blocked_range3d<size_t>(si,ei,sj,ej,sk,ek), 0.f,
                [&](const tbb::blocked_range3d<size_t> &r, float init)->T {
                    for(size_t i=r.pages().begin(), i_end=r.pages().end(); i<i_end; i++){
                        for (size_t j=r.rows().begin(), j_end=r.rows().end(); j<j_end; j++){
                            for (size_t k=r.cols().begin(), k_end=r.cols().end(); k<k_end; k++){
                                typedef std::vector<FluidParticle*>::iterator range_iterator;
                                init += tbb::parallel_reduce(
                                        tbb::blocked_range<range_iterator>(_MAC(i,j,k).begin(), _MAC(i,j,k).end()), 0.f,
                                        [&](const tbb::blocked_range<range_iterator> &r2, float init2)->T {
                                            for (range_iterator p = r2.begin(); p != r2.end(); p++) {
                                                FluidParticle* particle = *p;
                                                float dist = glm::distance(particle->pos, gridPos);
                                                float weight = kernel(dist, radius/2.f);
                                                T temp;
                                                void *address = (void *) particle + offset;
                                                std::memcpy(&temp, address, attributeSize);
                                                init2 += temp * (weight / totalWeight);
                                            }
                                            return init2;
                                        }, std::plus<T>()
                                );
                            }
                        }
                    }
                    return init;
                },
                std::plus<float>()
        );

        grid(I,J,K) = gridVal;

#else
        for (size_t i = si; i < ei; i++) {
            for (size_t j = sj; j < ej; j++) {
                for (size_t k = sk; k < ek; k++) {
                    for (FluidParticle const *particle : _MAC(i, j, k)) {
                        float dist = glm::distance(particle->pos, gridPos);
                        float weight = kernel(dist, radius / 2.f);
                        totalWeight += weight;
                    }
                }
            }
        }

        if (totalWeight == 0) {
            grid(I,J,K) = zeroVal;
            return;
        }

        T temp;
        T gridVal = zeroVal;
        for (size_t i = si; i < ei; i++) {
            for (size_t j = sj; j < ej; j++) {
                for (size_t k = sk; k < ek; k++) {
                    for (FluidParticle const *particle : _MAC(i, j, k)) {
                        float dist = glm::distance(particle->pos, gridPos);
                        float weight = kernel(dist, radius / 2.f);
                        void *address = (void *) particle + offset;
                        std::memcpy(&temp, address, attributeSize);
                        gridVal += temp * (weight / totalWeight);
                    }
                }
            }
        }

        grid(I,J,K) = gridVal;
#endif
    });

#endif // endif splatting
}

template<typename T> T FluidSolver::interpolateAttribute(const glm::vec3 &pos, Grid<T> &grid) {
    glm::vec3 idx = grid.fractionalIndexOf(pos);
    size_t i = (size_t) floor(idx.x);
    size_t j = (size_t) floor(idx.y);
    size_t k = (size_t) floor(idx.z);
    size_t I = (size_t) MATHIFELSE(ceil(idx.x), grid.countX()-1, ceil(idx.x) >= grid.countX());
    size_t J = (size_t) MATHIFELSE(ceil(idx.y), grid.countY()-1, ceil(idx.y) >= grid.countY());
    size_t K = (size_t) MATHIFELSE(ceil(idx.z), grid.countZ()-1, ceil(idx.z) >= grid.countZ());

    T k1, k2, k3, k4, j1, j2, val;

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
    enforceBoundary();
    gravitySolve(step);
    extrapolateVelocity();
    transferVelocitiesToParticles();
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