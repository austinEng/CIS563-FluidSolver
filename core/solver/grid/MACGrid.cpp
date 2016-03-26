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
        Grid<T>(origin, size*glm::vec3(0.5f,0.5f,0.5f), dim, size),
        _gU(Grid<float>(origin, size*glm::vec3(0.0f,0.5f,0.5f), dim, size)),
        _gV(Grid<float>(origin, size*glm::vec3(0.5f,0.0f,0.5f), dim, size)),
        _gW(Grid<float>(origin, size*glm::vec3(0.5f,0.5f,0.0f), dim, size)),
        _gU_old(Grid<float>(origin, size*glm::vec3(0.0f,0.5f,0.5f), dim, size)),
        _gV_old(Grid<float>(origin, size*glm::vec3(0.5f,0.0f,0.5f), dim, size)),
        _gW_old(Grid<float>(origin, size*glm::vec3(0.5f,0.5f,0.0f), dim, size)),
        _gP(Grid<float>(origin, size*glm::vec3(0.5f,0.5f,0.5f), dim, size)),
        _gType(Grid<int>(origin, size*glm::vec3(0.5f,0.5f,0.5f), dim, size)) {

}

template class MACGrid<std::vector<FluidParticle*, std::allocator<FluidParticle*> > >;