//
// Created by austin on 3/20/16.
//

#include "MACGrid.h"
#include <core/solver/FluidParticle.h>

template <typename T> MACGrid<T>::MACGrid() {

}

template <typename T> MACGrid<T>::~MACGrid() {

}

template <typename T> MACGrid<T>::MACGrid(const glm::vec3 &origin, const glm::vec3 &dim, float size) :
        Grid<T>(origin, dim, size),
        _gU(Grid<float>(size*glm::vec3(0.0f,0.5f,0.5f) + origin, dim, size)),
        _gV(Grid<float>(size*glm::vec3(0.5f,0.0f,0.5f) + origin, dim, size)),
        _gW(Grid<float>(size*glm::vec3(0.5f,0.5f,0.0f) + origin, dim, size)),
        _gU_old(Grid<float>(size*glm::vec3(0.0f,0.5f,0.5f) + origin, dim, size)),
        _gV_old(Grid<float>(size*glm::vec3(0.5f,0.0f,0.5f) + origin, dim, size)),
        _gW_old(Grid<float>(size*glm::vec3(0.5f,0.5f,0.0f) + origin, dim, size)),
        _gP(Grid<float>(size*glm::vec3(0.5f,0.5f,0.5f) + origin, dim, size)) {

}

template class MACGrid<std::vector<FluidParticle*, std::allocator<FluidParticle*> > >;