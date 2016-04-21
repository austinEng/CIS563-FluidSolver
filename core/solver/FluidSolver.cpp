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
#include <Eigen/Sparse>
#include <iostream>


float FluidSolver::g = -9.80665f;

FluidSolver::FluidSolver(float particleSep, float gridSize) : particle_radius(particleSep), _cell_size(gridSize), frame(0) {
//    particle_radius /= 2;
}

FluidSolver::~FluidSolver() {
    delete _container;
    delete _MAC;
}

void FluidSolver::setContainer(GeoObject* container) {
    _container = container;
    glm::vec3 size = _container->bound().dim();
    glm::vec3 origin = (_container->bound().center() - size / 2.f);
    _MAC = new MACGrid<std::vector<FluidParticle*> >(
            origin - glm::vec3(_cell_size, _cell_size, _cell_size),
            size + 2.f*glm::vec3(_cell_size, _cell_size, _cell_size),
            _cell_size
    );
//    _MAC = new MACGrid<std::vector<FluidParticle*> >(origin, size, _cell_size);

    std::function<void(size_t, size_t, size_t)> setSolid = [&](size_t i, size_t j, size_t k) {
        _MAC->_gType(i,j,k) = SOLID;
    };

    _MAC->_gType.iterateRegion(0,0,0, 1,_MAC->_gType.countY(),_MAC->_gType.countZ(), setSolid);
    _MAC->_gType.iterateRegion(_MAC->_gType.countX()-1,0,0, _MAC->_gType.countX(),_MAC->_gType.countY(),_MAC->_gType.countZ(), setSolid);
    _MAC->_gType.iterateRegion(0,0,0, _MAC->_gType.countX(),1,_MAC->_gType.countZ(), setSolid);
    _MAC->_gType.iterateRegion(0,_MAC->_gType.countY()-1,0, _MAC->_gType.countX(),_MAC->_gType.countY(),_MAC->_gType.countZ(), setSolid);
    _MAC->_gType.iterateRegion(0,0,0, _MAC->_gType.countX(),_MAC->_gType.countY(),1, setSolid);
    _MAC->_gType.iterateRegion(0,0,_MAC->_gType.countZ()-1, _MAC->_gType.countX(),_MAC->_gType.countY(),_MAC->_gType.countZ(), setSolid);

    /*_MAC = MACGrid<std::vector<FluidParticle*> >(
            origin,
            size,
            _cell_size
    );*/
}

/*
 * Loop over fluid bounds to generate particles
 */
void FluidSolver::addFluid(const GeoObject &fluid) {
//    FluidParticle p;
//
//    p.pos = glm::vec3(-1,5,1);
//    _particles.push_back(p);
//    _MAC->_gType.at(p.pos) = FLUID;
//
//    p.pos = glm::vec3(-1,3,1);
//    _particles.push_back(p);
//    _MAC->_gType.at(p.pos) = FLUID;
//
//    p.pos = glm::vec3(-1,5,3);
//    _particles.push_back(p);
//    _MAC->_gType.at(p.pos) = FLUID;
//
//    p.pos = glm::vec3(-1,3,3);
//    _particles.push_back(p);
//    _MAC->_gType.at(p.pos) = FLUID;
//
//    p.pos = glm::vec3(-3,5,1);
//    _particles.push_back(p);
//    _MAC->_gType.at(p.pos) = FLUID;
//
//    p.pos = glm::vec3(-3,3,1);
//    _particles.push_back(p);
//    _MAC->_gType.at(p.pos) = FLUID;
//
//    p.pos = glm::vec3(-3,5,3);
//    _particles.push_back(p);
//    _MAC->_gType.at(p.pos) = FLUID;
//
//    p.pos = glm::vec3(-3,3,3);
//    _particles.push_back(p);
//    _MAC->_gType.at(p.pos) = FLUID;
//    return;

    const Bound& b = fluid.bound();
    for (float x = b.minX(); x < b.maxX(); x += particle_radius) {
        for (float y = b.minY(); y < b.maxY(); y += particle_radius) {
            for (float z = b.minZ(); z < b.maxZ(); z += particle_radius) {
                glm::vec3 pos = glm::vec3(x, y, z) + glm::vec3(_cell_size)/2.f;
                if (fluid.contains(pos)) {
                    FluidParticle p;
                    p.pos = pos;
                    _particles.push_back(p);
                    _MAC->_gType.at(pos) = FLUID;
                }
            }
        }
    }
    std::cout << "Added " << _particles.size() << " particles" << std::endl;
}

void FluidSolver::init() {

    for (FluidParticle &particle : _particles) {
        particle.cell = _MAC->indexOf(particle.pos);
        _MAC->atIdx(particle.cell).push_back(&particle);
    }
    _MAC->_gU.clear(0);
    _MAC->_gV.clear(0);
    _MAC->_gW.clear(0);

//    for (size_t idx = 0; idx < _MAC->_gU.size(); idx++) {
//        size_t i, j, k;
//        _MAC->_gU.toIJK(idx, i,j,k);
//        std::cout << idx << "; " << i << "," << j << "," << k << "; " << _MAC->_gU.fromIJK(i,j,k) << std::endl;
//    }
//
//    std::size_t velOffset = offsetof(FluidParticle, vel);
//    std::size_t U_offset = velOffset + offsetof(glm::vec3, x);
//    std::size_t V_offset = velOffset + offsetof(glm::vec3, y);
//    std::size_t W_offset = velOffset + offsetof(glm::vec3, z);
//    particleAttributeToGrid(U_offset, _MAC->_gU_old, _cell_size, 0.f);
//    particleAttributeToGrid(V_offset, _MAC->_gV_old, _cell_size, 0.f);
//    particleAttributeToGrid(W_offset, _MAC->_gW_old, _cell_size, 0.f);
//    particleAttributeToGrid(U_offset, _MAC->_gU, _cell_size, 0.f);
//    particleAttributeToGrid(V_offset, _MAC->_gV, _cell_size, 0.f);
//    particleAttributeToGrid(W_offset, _MAC->_gW, _cell_size, 0.f);
}

void FluidSolver::projectVelocitiesToGrid() {
    std::size_t velOffset = offsetof(FluidParticle, vel);
    std::size_t U_offset = velOffset + offsetof(glm::vec3, x);
    std::size_t V_offset = velOffset + offsetof(glm::vec3, y);
    std::size_t W_offset = velOffset + offsetof(glm::vec3, z);

#ifdef USETBB
    tbb::parallel_invoke(
            [&](){
                particleAttributeToGrid(U_offset, _MAC->_gU, _cell_size, 0.f);
            },
            [&](){
                particleAttributeToGrid(V_offset, _MAC->_gV, _cell_size, 0.f);
            },
            [&](){
                particleAttributeToGrid(W_offset, _MAC->_gW, _cell_size, 0.f);
            }
    );
    tbb::parallel_invoke(
            [&]() {
                _MAC->_gU_old = _MAC->_gU;
            },
            [&]() {
                _MAC->_gV_old = _MAC->_gV;
            },
            [&]() {
                _MAC->_gW_old = _MAC->_gW;
            }
    );

#else
    particleAttributeToGrid(U_offset, _MAC->_gU, _cell_size, 0.f);
    particleAttributeToGrid(V_offset, _MAC->_gV, _cell_size, 0.f);
    particleAttributeToGrid(W_offset, _MAC->_gW, _cell_size, 0.f);

    _MAC->_gU_old = _MAC->_gU;
    _MAC->_gV_old = _MAC->_gV;
    _MAC->_gW_old = _MAC->_gW;
#endif
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
                                          float vel = interpolateAttribute(particle.pos, _MAC->_gU);
                                          float oldVel = interpolateAttribute(particle.pos, _MAC->_gU_old);
                                          particle.vel.x = vel*smooth + (particle.vel.x +(vel - oldVel))*(1.f-smooth);
                                      }
                                  });
            },
            [&]() {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, _particles.size()),
                                  [&](const tbb::blocked_range<size_t> &r) {
                                      for (size_t i = r.begin(); i != r.end(); ++i) {
                                          FluidParticle &particle = _particles[i];
                                          float vel = interpolateAttribute(particle.pos, _MAC->_gV);
                                          float oldVel = interpolateAttribute(particle.pos, _MAC->_gV_old);
                                          particle.vel.y = vel*smooth + (particle.vel.y +(vel - oldVel))*(1.f-smooth);
                                      }
                                  });
            },
            [&]() {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, _particles.size()),
                                  [&](const tbb::blocked_range<size_t> &r) {
                                      for (size_t i = r.begin(); i != r.end(); ++i) {
                                          FluidParticle &particle = _particles[i];
                                          float vel = interpolateAttribute(particle.pos, _MAC->_gW);
                                          float oldVel = interpolateAttribute(particle.pos, _MAC->_gW_old);
                                          particle.vel.z = vel*smooth + (particle.vel.z +(vel - oldVel))*(1.f-smooth);
                                      }
                                  });
            }
    );
#else
    for (FluidParticle &particle : _particles) {
        float vel = interpolateAttribute(particle.pos, _MAC->_gU);
        float oldVel = interpolateAttribute(particle.pos, _MAC->_gU_old);
        particle.vel.x = vel*smooth + (particle.vel.x +(vel - oldVel))*(1.f-smooth);
    }
    for (FluidParticle &particle : _particles) {
        float vel = interpolateAttribute(particle.pos, _MAC->_gV);
        float oldVel = interpolateAttribute(particle.pos, _MAC->_gV_old);
        particle.vel.y = vel*smooth + (particle.vel.y +(vel - oldVel))*(1.f-smooth);
    }
    for (FluidParticle &particle : _particles) {
        float vel = interpolateAttribute(particle.pos, _MAC->_gW);
        float oldVel = interpolateAttribute(particle.pos, _MAC->_gW_old);
        particle.vel.z = vel*smooth + (particle.vel.z +(vel - oldVel))*(1.f-smooth);
    }
#endif
}

void FluidSolver::enforceBoundary() {
    _MAC->_gType.iterate([&](size_t i, size_t j, size_t k) {
        switch (_MAC->_gType(i,j,k)) {
            case EMPTY:break;
            case FLUID:break;
            case SOLID:
                /*_MAC->_gU(i,j,k) = 0;
                _MAC->_gV(i,j,k) = 0;
                _MAC->_gW(i,j,k) = 0;
                _MAC->_gU(i+1,j,k) = 0;
                _MAC->_gV(i,j+1,k) = 0;
                _MAC->_gW(i,j,k+1) = 0;
                */
                if (i == 0 || _MAC->_gType(i-1,j,k) != SOLID) {
                    _MAC->_gU(i,j,k) = std::min(0.f, _MAC->_gU(i,j,k));
                }
                if (i == _MAC->_gType.countX() - 1 || _MAC->_gType(i+1,j,k) != SOLID) {
                    _MAC->_gU(i+1,j,k) = std::max(0.f, _MAC->_gU(i+1,j,k));
                }
                if (j == 0 || _MAC->_gType(i,j-1,k) != SOLID) {
                    _MAC->_gV(i,j,k) = std::min(0.f, _MAC->_gV(i,j,k));
                }
                if (j == _MAC->_gType.countY() - 1 || _MAC->_gType(i,j+1,k) != SOLID) {
                    _MAC->_gV(i,j+1,k) = std::max(0.f, _MAC->_gV(i,j+1,k));
                }
                if (k == 0 || _MAC->_gType(i,j,k-1) != SOLID) {
                    _MAC->_gW(i,j,k) = std::min(0.f, _MAC->_gW(i,j,k));
                }
                if (k == _MAC->_gType.countZ() - 1 || _MAC->_gType(i,j,k+1) != SOLID) {
                    _MAC->_gW(i,j,k+1) = std::max(0.f, _MAC->_gW(i,j,k+1));
                }
                break;
            default:break;
        }
    });
    //return;


    /*_MAC->_gU.iterateRegion(0,0,0, 1,_MAC->_gU.countY(),_MAC->_gU.countZ(), [&](size_t i, size_t j, size_t k) {
        _MAC->_gU(i,j,k) = std::min(0.f, _MAC->_gU(i,j,k));
    });
    _MAC->_gU.iterateRegion(_MAC->_gU.countX()-1,0,0, _MAC->_gU.countX(),_MAC->_gU.countY(),_MAC->_gU.countZ(), [&](size_t i, size_t j, size_t k) {
        _MAC->_gU(i,j,k) = std::max(0.f, _MAC->_gU(i,j,k));
    });
    _MAC->_gV.iterateRegion(0,0,0, _MAC->_gV.countX(),1,_MAC->_gV.countZ(), [&](size_t i, size_t j, size_t k) {
        _MAC->_gV(i,j,k) = 0;
    });
    _MAC->_gV.iterateRegion(0,_MAC->_gV.countY()-1,0, _MAC->_gV.countX(),_MAC->_gV.countY(),_MAC->_gV.countZ(), [&](size_t i, size_t j, size_t k) {
        _MAC->_gV(i,j,k) = 0;
    });
    _MAC->_gW.iterateRegion(0,0,0, _MAC->_gW.countX(),_MAC->_gW.countY(),1, [&](size_t i, size_t j, size_t k) {
        _MAC->_gW(i,j,k) = 0;
    });
    _MAC->_gW.iterateRegion(0,0,_MAC->_gW.countZ()-1, _MAC->_gW.countX(),_MAC->_gW.countY(),_MAC->_gW.countZ(), [&](size_t i, size_t j, size_t k) {
        _MAC->_gW(i,j,k) = 0;
    });*/
}

inline void pressureMatrixHelper(std::vector<Eigen::Triplet<double>> &coeffs, const Grid<int> &grid, size_t &IDX,
                                 int &count, const float &scale, size_t i, size_t j, size_t k) {
    size_t idx = grid.fromIJK(i,j,k);
    if (grid(i,j,k) == FLUID || grid(i,j,k) == EMPTY) {
        count++;
        if (grid(i,j,k) == FLUID) {
            coeffs.push_back(Eigen::Triplet<double>(IDX, idx, -scale));
        }
    }
}

void FluidSolver::pressureSolve(float step) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> coefficientsA;
    std::vector<T> coefficientsB;
    Eigen::SparseMatrix<double> A(_MAC->_gType.size(), _MAC->_gType.size());
    Eigen::SparseMatrix<double> b(_MAC->_gType.size(), 1);
    Eigen::SparseVector<double> x(_MAC->_gType.size());
    A.setZero();
    b.setZero();
    x.setZero();

    float scale = step / (1.f*_cell_size*_cell_size);

    _MAC->_gType.iterate([&](size_t I, size_t J, size_t K) {
        size_t IDX = _MAC->_gType.fromIJK(I,J,K);
        if (_MAC->_gType(I,J,K) == FLUID) {
            int count = 0;

            if (I > 0) {                            // if I-1 >= 0
                pressureMatrixHelper(coefficientsA, _MAC->_gType, IDX, count, scale, I-1,J,K);
            }
            if (I + 1 < _MAC->_gType.countX()) {     // if I + 1 < countX
                pressureMatrixHelper(coefficientsA, _MAC->_gType, IDX, count, scale, I+1,J,K);
            }
            if (J > 0) {                            // if J-1 >= 0
                pressureMatrixHelper(coefficientsA, _MAC->_gType, IDX, count, scale, I,J-1,K);
            }
            if (J + 1 < _MAC->_gType.countY()) {     // if J + 1 < countY
                pressureMatrixHelper(coefficientsA, _MAC->_gType, IDX, count, scale, I,J+1,K);
            }
            if (K > 0) {                            // if K-1 >= 0
                pressureMatrixHelper(coefficientsA, _MAC->_gType, IDX, count, scale, I,J,K-1);
            }
            if (K + 1 < _MAC->_gType.countZ()) {     // if K + 1 < countZ
                pressureMatrixHelper(coefficientsA, _MAC->_gType, IDX, count, scale, I,J,K+1);
            }

            coefficientsA.push_back(T(IDX, IDX, count*scale));

            float div =
            -(_MAC->_gU(I+1,J,K) - _MAC->_gU(I,J,K)) / _cell_size +
            -(_MAC->_gV(I,J+1,K) - _MAC->_gV(I,J,K)) / _cell_size +
            -(_MAC->_gW(I,J,K+1) - _MAC->_gW(I,J,K)) / _cell_size;

            if (I == 0 || _MAC->_gType(I-1,J,K) == SOLID) {
                div -= (_MAC->_gU(I,J,K) - 0) / _cell_size;
            }
            if (I == _MAC->_gType.countX() - 1 || _MAC->_gType(I+1,J,K) == SOLID) {
                div -= (0 - _MAC->_gU(I+1,J,K)) / _cell_size;
            }
            if (J == 0 || _MAC->_gType(I,J-1,K) == SOLID) {
                div -= (_MAC->_gV(I,J,K) - 0) / _cell_size;
            }
            if (J == _MAC->_gType.countY() - 1 || _MAC->_gType(I,J+1,K) == SOLID) {
                div -= (0 - _MAC->_gV(I,J+1,K)) / _cell_size;
            }
            if (K == 0 || _MAC->_gType(I,J,K-1) == SOLID) {
                div -= (_MAC->_gW(I,J,K) - 0) / _cell_size;
            }
            if (K == _MAC->_gType.countZ() - 1 || _MAC->_gType(I,J,K+1) == SOLID) {
                div -= (0 - _MAC->_gW(I,J,K+1)) / _cell_size;
            }

            coefficientsB.push_back(T(IDX,0,div));
        }
    }, false);

    A.setFromTriplets(coefficientsA.begin(), coefficientsA.end());

    b.setFromTriplets(coefficientsB.begin(), coefficientsB.end());
//    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::IncompleteCholesky<double> > cg(A);
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::IdentityPreconditioner> cg(A);
    x = cg.solve(b);

    _MAC->_gP.clear(0);
    for (Eigen::SparseVector<double>::InnerIterator it(x); it; ++it) {
        if (_MAC->_gType(it.index()) == FLUID) {
            _MAC->_gP(it.index()) = it.value();
        }
    }

    scale = step/(1.f*_cell_size);
    _MAC->_gU.iterate([&](size_t i, size_t j, size_t k) {
        bool leftExists = i > 0;
        bool rightExists = i < _MAC->_gP.countX();
        bool leftFluid = leftExists && _MAC->_gType(i-1,j,k) == FLUID;
        bool rightFluid = rightExists && _MAC->_gType(i,j,k) == FLUID;
        if ((leftFluid || rightFluid)) {
            float delP = _MAC->_gP(i,j,k) - _MAC->_gP(i-1,j,k);
            _MAC->_gU(i,j,k) -= scale * delP;
        }
    });

    _MAC->_gV.iterate([&](size_t i, size_t j, size_t k) {
        bool leftExists = j > 0;
        bool rightExists = j < _MAC->_gP.countY();
        bool leftFluid = leftExists && _MAC->_gType(i,j-1,k) == FLUID;
        bool rightFluid = rightExists && _MAC->_gType(i,j,k) == FLUID;
        if ((leftFluid || rightFluid)) {
            float delP = _MAC->_gP(i,j,k) - _MAC->_gP(i,j-1,k);
            _MAC->_gV(i,j,k) -= scale * delP;
        }
    });

    _MAC->_gW.iterate([&](size_t i, size_t j, size_t k) {
        bool leftExists = k > 0;
        bool rightExists = k < _MAC->_gP.countZ();
        bool leftFluid = leftExists && _MAC->_gType(i,j,k-1) == FLUID;
        bool rightFluid = rightExists && _MAC->_gType(i,j,k) == FLUID;
        if ((leftFluid || rightFluid)) {
            float delP = _MAC->_gP(i,j,k) - _MAC->_gP(i,j,k-1);
            _MAC->_gW(i,j,k) -= scale * delP;
        }
    });

}

void FluidSolver::gravitySolve(float step) {
    _MAC->_gV.iterate([&](size_t i, size_t j, size_t k) {
        _MAC->_gV(i,j,k) += g*step;
    });
}

void FluidSolver::extrapolateVelocity() {
    _MAC->_gU.iterate([&](size_t i, size_t j, size_t k) {
        bool shouldExtrapolate = (_MAC->_gType.checkIdx(i-1,j,k) && _MAC->_gType(i-1,j,k) != FLUID) || (_MAC->_gType.checkIdx(i,j,k) && _MAC->_gType(i,j,k) != FLUID);
        if (shouldExtrapolate) {
            bool fromUp = (_MAC->_gType.checkIdx(i-1,j+1,k) && _MAC->_gType(i-1,j+1,k) == FLUID) || (_MAC->_gType.checkIdx(i,j+1,k) && _MAC->_gType(i,j+1,k) == FLUID);
            bool fromDown = (_MAC->_gType.checkIdx(i-1,j-1,k) && _MAC->_gType(i-1,j-1,k) == FLUID) || (_MAC->_gType.checkIdx(i,j-1,k) && _MAC->_gType(i,j-1,k) == FLUID);
            bool fromFront = (_MAC->_gType.checkIdx(i-1,j,k+1) && _MAC->_gType(i-1,j,k+1) == FLUID) || (_MAC->_gType.checkIdx(i,j,k+1) && _MAC->_gType(i,j,k+1) == FLUID);
            bool fromBack = (_MAC->_gType.checkIdx(i-1,j,k-1) && _MAC->_gType(i-1,j,k-1) == FLUID) || (_MAC->_gType.checkIdx(i,j,k-1) && _MAC->_gType(i,j,k-1) == FLUID);

            float val = 0;
            int count = fromUp + fromDown + fromFront + fromBack;
            if (fromUp) val += _MAC->_gU(i,j+1,k);
            if (fromDown) val += _MAC->_gU(i,j-1,k);
            if (fromFront) val += _MAC->_gU(i,j,k+1);
            if (fromBack) val += _MAC->_gU(i,j,k-1);

            if (count > 0) {
                _MAC->_gU(i,j,k) = val / count;
            }
        }
    });

    _MAC->_gV.iterate([&](size_t i, size_t j, size_t k) {
        bool shouldExtrapolate = (_MAC->_gType.checkIdx(i,j-1,k) && _MAC->_gType(i,j-1,k) != FLUID) || (_MAC->_gType.checkIdx(i,j,k) && _MAC->_gType(i,j,k) != FLUID);
        if (shouldExtrapolate) {
            bool fromUp = (_MAC->_gType.checkIdx(i+1,j-1,k) && _MAC->_gType(i+1,j-1,k) == FLUID) || (_MAC->_gType.checkIdx(i+1,j,k) && _MAC->_gType(i+1,j,k) == FLUID);
            bool fromDown = (_MAC->_gType.checkIdx(i-1,j-1,k) && _MAC->_gType(i-1,j-1,k) == FLUID) || (_MAC->_gType.checkIdx(i-1,j,k) && _MAC->_gType(i-1,j,k) == FLUID);
            bool fromFront = (_MAC->_gType.checkIdx(i,j-1,k+1) && _MAC->_gType(i,j-1,k+1) == FLUID) || (_MAC->_gType.checkIdx(i,j,k+1) && _MAC->_gType(i,j,k+1) == FLUID);
            bool fromBack = (_MAC->_gType.checkIdx(i,j-1,k-1) && _MAC->_gType(i,j-1,k-1) == FLUID) || (_MAC->_gType.checkIdx(i,j,k-1) && _MAC->_gType(i,j,k-1) == FLUID);

            float val = 0;
            int count = fromUp + fromDown + fromFront + fromBack;
            if (fromUp) val += _MAC->_gV(i+1,j,k);
            if (fromDown) val += _MAC->_gV(i-1,j,k);
            if (fromFront) val += _MAC->_gV(i,j,k+1);
            if (fromBack) val += _MAC->_gV(i,j,k-1);

            if (count > 0) {
                _MAC->_gV(i,j,k) = val / count;
            }
        }
    });

    _MAC->_gW.iterate([&](size_t i, size_t j, size_t k) {
        bool shouldExtrapolate = (_MAC->_gType.checkIdx(i,j,k-1) && _MAC->_gType(i,j,k-1) != FLUID) || (_MAC->_gType.checkIdx(i,j,k) && _MAC->_gType(i,j,k) != FLUID);
        if (shouldExtrapolate) {
            bool fromUp = (_MAC->_gType.checkIdx(i+1,j,k-1) && _MAC->_gType(i+1,j,k-1) == FLUID) || (_MAC->_gType.checkIdx(i+1,j,k) && _MAC->_gType(i+1,j,k) == FLUID);
            bool fromDown = (_MAC->_gType.checkIdx(i-1,j,k-1) && _MAC->_gType(i-1,j,k-1) == FLUID) || (_MAC->_gType.checkIdx(i-1,j,k) && _MAC->_gType(i-1,j,k) == FLUID);
            bool fromFront = (_MAC->_gType.checkIdx(i,j+1,k-1) && _MAC->_gType(i,j+1,k-1) == FLUID) || (_MAC->_gType.checkIdx(i,j+1,k) && _MAC->_gType(i,j+1,k) == FLUID);
            bool fromBack = (_MAC->_gType.checkIdx(i,j-1,k-1) && _MAC->_gType(i,j-1,k-1) == FLUID) || (_MAC->_gType.checkIdx(i,j-1,k) && _MAC->_gType(i,j-1,k) == FLUID);

            float val = 0;
            int count = fromUp + fromDown + fromFront + fromBack;
            if (fromUp) val += _MAC->_gW(i+1,j,k);
            if (fromDown) val += _MAC->_gW(i-1,j,k);
            if (fromFront) val += _MAC->_gW(i,j+1,k);
            if (fromBack) val += _MAC->_gW(i,j-1,k);

            if (count > 0) {
                _MAC->_gW(i,j,k) = val / count;
            }
        }
    });

}

void FluidSolver::updateParticlePositions(float step) {

#ifdef USETBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, _particles.size()), [&](const tbb::blocked_range<size_t> &r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            FluidParticle &particle = _particles[i];
            particle.pos_old = particle.pos;

            glm::vec3 k1 = step*particle.vel;
            particle.pos += step * glm::vec3(
                    interpolateAttribute(particle.pos + 0.5f*k1, _MAC->_gU),
                    interpolateAttribute(particle.pos + 0.5f*k1, _MAC->_gV),
                    interpolateAttribute(particle.pos + 0.5f*k1, _MAC->_gW)
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
                //particle.col = glm::vec3(1,0,0);
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
            //particle.vel *= mask;
            particle.pos = particle.pos_old;
        }
    }
#endif
}

void FluidSolver::updateCells() {
    _MAC->clear(std::vector<FluidParticle*>());
    _MAC->_gType.clear(EMPTY);
    for (FluidParticle &particle : _particles) {
        particle.cell = _MAC->indexOf(particle.pos);
        _MAC->_gType(particle.cell) = FLUID;
        if (_MAC->checkIdx(particle.cell)) {
            _MAC->atIdx(particle.cell).push_back(&particle);
        } else {
            //std::cerr << "particle out of bounds" << std::endl;
        }

    }

    std::function<void(size_t, size_t, size_t)> setSolid = [&](size_t i, size_t j, size_t k) {
        _MAC->_gType(i,j,k) = SOLID;
    };

    _MAC->_gType.iterateRegion(0,0,0, 1,_MAC->_gType.countY(),_MAC->_gType.countZ(), setSolid);
    _MAC->_gType.iterateRegion(_MAC->_gType.countX()-1,0,0, _MAC->_gType.countX(),_MAC->_gType.countY(),_MAC->_gType.countZ(), setSolid);
    _MAC->_gType.iterateRegion(0,0,0, _MAC->_gType.countX(),1,_MAC->_gType.countZ(), setSolid);
    _MAC->_gType.iterateRegion(0,_MAC->_gType.countY()-1,0, _MAC->_gType.countX(),_MAC->_gType.countY(),_MAC->_gType.countZ(), setSolid);
    _MAC->_gType.iterateRegion(0,0,0, _MAC->_gType.countX(),_MAC->_gType.countY(),1, setSolid);
    _MAC->_gType.iterateRegion(0,0,_MAC->_gType.countZ()-1, _MAC->_gType.countX(),_MAC->_gType.countY(),_MAC->_gType.countZ(), setSolid);

}

inline float kernel(float r, float h) {
//    float e = 2*r/h;
//    return 1.f/(PI*h*h*h) * MATHIFELSE(
//            1.f - 3.f/2.f * e*e  + 3.f/4.f * e*e*e,
//            MATHIFELSE(
//                    1.f/4.f * (2-e)*(2-e)*(2-e),
//                    0,
//                    e > 2
//            ),
//            e > 1
//    );
//    return (float) (315 / (64 * PI * pow(h, 9)) * MATHIFELSE(0, pow(h * h - r * r, 3), 0 <= r && r <= h));
    return MATHIFELSE(0, 1.f-r/h, 0 <= r && r <= h);
}


template<typename T> void FluidSolver::particleAttributeToGrid(std::size_t offset, Grid<T> &grid, float radius, T zeroVal) {
    std::size_t attributeSize = sizeof(T);
    std::size_t cellRadius = (size_t) glm::ceil(radius / _cell_size);

//    grid.clear(zeroVal);
//    std::vector<float> weights(grid.countX() * grid.countY() * grid.countZ());
//
//    iterParticles([&](FluidParticle &particle) {
//        size_t I,J,K;
//        grid.indexOf(particle.pos, I, J, K);
//        glm::vec3 gridPos = grid.positionOf(I,J,K);
//
//        grid.iterateNeighborhood(I,J,K,cellRadius, [&](size_t i, size_t j, size_t k) {
//            float dist = glm::distance(particle.pos, gridPos);
//            size_t idx = grid.fromIJK(i,j,k);
//            weights[idx] += kernel(dist, radius);
//        });
//    }, false);
//
//    iterParticles([&](FluidParticle &particle) {
//        size_t I,J,K;
//        grid.indexOf(particle.pos, I, J, K);
//        glm::vec3 gridPos = grid.positionOf(I,J,K);
//
//        grid.iterateNeighborhood(I,J,K,cellRadius, [&](size_t i, size_t j, size_t k) {
//            float dist = glm::distance(particle.pos, gridPos);
//            size_t idx = grid.fromIJK(i,j,k);
//            T temp;
//            void *address = (void *) &particle + offset;
//            std::memcpy(&temp, address, attributeSize);
//            grid(i,j,k) += temp * (kernel(dist, radius) / weights[idx]);
//        });
//    }, false);
//
//    return;

    grid.iterate([&](size_t I, size_t J, size_t K) {
        glm::vec3 gridPos = grid.positionOf(I,J,K);

        size_t mI, mJ, mK, si, ei, sj, ej, sk, ek;
        _MAC->indexOf(gridPos, mI, mJ, mK);
        _MAC->getNeighboorhood(mI, mJ, mK, cellRadius, si, ei, sj, ej, sk, ek);

        float totalWeight = 0.f;

//        totalWeight = tbb::parallel_reduce(tbb::blocked_range3d<size_t>(si, ei, sj, ej, sk, ek), 0.f, [&](const tbb::blocked_range3d<size_t> &r, float init)->float {
//            for (size_t i = r.rows().begin(); i < r.rows().end(); i++) {
//                for (size_t j = r.cols().begin(); j < r.cols().end(); j++) {
//                    for (size_t k = r.pages().begin(); k < r.pages().end(); k++) {
//                        for (FluidParticle const *particle : _MAC(i, j, k)) {
//                            float dist = glm::distance(particle->pos, gridPos);
//                            float weight = kernel(dist, 2*radius);
//                            init += weight;
//                        }
//                    }
//                }
//            }
//            return init;
//        }, std::plus<float>()
//        );

        for (size_t i = si; i < ei; i++) {
            for (size_t j = sj; j < ej; j++) {
                for (size_t k = sk; k < ek; k++) {
                    for (FluidParticle const *particle : _MAC->atIdx(i, j, k)) {
                        float dist = glm::distance2(particle->pos, gridPos);
                        float weight = kernel(dist, 2*radius*radius);
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
                    for (FluidParticle const *particle : _MAC->atIdx(i, j, k)) {
                        float dist = glm::distance2(particle->pos, gridPos);
                        float weight = kernel(dist, 2*radius*radius);
                        void *address = (void *) particle + offset;
                        std::memcpy(&temp, address, attributeSize);
                        gridVal += temp * (weight / totalWeight);
                    }
                }
            }
        }

        grid(I,J,K) = gridVal;
    });

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
    gravitySolve(step);
    enforceBoundary();
    pressureSolve(step);
    enforceBoundary();
    extrapolateVelocity();
    enforceBoundary();
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