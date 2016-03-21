//
// Created by austin on 3/20/16.
//

#include "Grid.h"
#include <core/solver/FluidParticle.h>

template <typename T> Grid<T>::Grid() {

}

template <typename T> Grid<T>::~Grid() {

}

template <typename T> Grid<T>::Grid(const glm::vec3 &offset, const glm::vec3 &dim, float size) :
        _offset(offset),
        _dim(dim),
        _cellSize(size),
        _cellCount(glm::ceil((_dim - _offset) / _cellSize)) {
    _contents = std::vector<T>((unsigned long) (_cellCount.x * _cellCount.y * _cellCount.z));
}

template <typename T> T& Grid<T>::operator()(std::size_t idx) {
    return _contents[idx];
}

template <typename T> const T& Grid<T>::operator()(std::size_t idx) const {
    return _contents[idx];
}

template <typename T> T& Grid<T>::operator()(std::size_t i, std::size_t j, std::size_t k) {
    return _contents[k*_cellCount.x*_cellCount.y + j*_cellCount.x + i];
}

template <typename T> const T& Grid<T>::operator()(std::size_t i, std::size_t j, std::size_t k) const {
    return _contents[k*_cellCount.x*_cellCount.y + j*_cellCount.x + i];
}


template <typename T> T& Grid<T>::atIdx(std::size_t i, std::size_t j, std::size_t k) {
    return _contents[k*_cellCount.x*_cellCount.y + j*_cellCount.x + i];
}

template <typename T> const T& Grid<T>::atIdx(std::size_t i, std::size_t j, std::size_t k) const {
    return _contents[k*_cellCount.x*_cellCount.y + j*_cellCount.x + i];
}


template <typename T> T& Grid<T>::operator()(const glm::ivec3 &idx) {
    return _contents[idx.z*_cellCount.x*_cellCount.y + idx.y*_cellCount.x + idx.x];
}

template <typename T> const T& Grid<T>::operator()(const glm::ivec3 &idx) const {
    return _contents[idx.z*_cellCount.x*_cellCount.y + idx.y*_cellCount.x + idx.x];;
}

template <typename T> T& Grid<T>::atIdx(const glm::ivec3 &idx) {
    return _contents[idx.z*_cellCount.x*_cellCount.y + idx.y*_cellCount.x + idx.x];
}

template <typename T> const T& Grid<T>::atIdx(const glm::ivec3 &idx) const {
    return _contents[idx.z*_cellCount.x*_cellCount.y + idx.y*_cellCount.x + idx.x];
}


template <typename T> T& Grid<T>::at(float x, float y, float z) {
    return at(glm::vec3(x, y, z));
}

template <typename T> const T& Grid<T>::at(float x, float y, float z) const {
    return at(glm::vec3(x, y, z));;
}

template <typename T> T& Grid<T>::at(const glm::vec3 &pos) {
    glm::vec3 indices = (pos - _offset) / _cellSize;
    return this->operator()((size_t) indices.x, (size_t) indices.y, (size_t) indices.z);
}

template <typename T> const T& Grid<T>::at(const glm::vec3 &pos) const {
    glm::vec3 indices = (pos - _offset) / _cellSize;
    return this->operator()((size_t) indices.x, (size_t) indices.y, (size_t) indices.z);
}

template <typename T> glm::ivec3 Grid<T>::indexOf(const glm::vec3 &pos) const {
    glm::vec3 indices = (pos - _offset) / _cellSize;
    int i = (int) indices.x;
    int j = (int) indices.y;
    int k = (int) indices.z;
    if (i >= _cellCount.x ) i = -1;
    if (j >= _cellCount.y ) j = -1;
    if (k >= _cellCount.z ) k = -1;
    return glm::ivec3(i, j, k);
}

template <typename T> glm::vec3 Grid<T>::positionOf(const glm::ivec3 &idx) const {
    return glm::vec3(idx.x * _cellSize, idx.y * _cellSize, idx.z * _cellSize) + _offset;
}

template <typename T> glm::vec3 Grid<T>::positionOf(size_t i, size_t j, size_t k) const {
    return glm::vec3(i * _cellSize, j * _cellSize, k * _cellSize) + _offset;
}

template <typename T> glm::vec3 Grid<T>::fractionalIndexOf(const glm::vec3 &pos) const {
    return (pos - _offset) / _cellSize;
}

template <typename T> glm::ivec3 Grid<T>::toIJK(const std::size_t index) const {
    int i = (int) (index % _cellCount.z);
    int j = (int) ((index / _cellCount.z) & _cellCount.y);
    int k = (int) (index / (_cellCount.y * _cellCount.z));
    return glm::ivec3(i,j,k);
}

template <typename T> std::size_t Grid<T>::fromIJK(const std::size_t i, const std::size_t j, const std::size_t k) const {
    return (size_t) (k * _cellCount.x * _cellCount.y + j * _cellCount.x + i);
}

template <typename T> std::size_t Grid<T>::fromIJK(const glm::ivec3 &ijk) const {
    return (size_t) (ijk.z * _cellCount.x * _cellCount.y + ijk.y * _cellCount.x + ijk.x);
}

template <typename T> void Grid<T>::iterate(const std::function<void(size_t i, size_t j, size_t k)> &cb) {
    for (size_t i = 0; i < _cellCount.x; i++) {
        for (size_t j = 0; j < _cellCount.y; j++) {
            for (size_t k = 0; k < _cellCount.z; k++) {
                cb(i,j,k);
            }
        }
    }
}

template <typename T> void Grid<T>::clear(const T &zeroVal) {
    for (size_t i = 0; i < _contents.size(); i++) {
        _contents[i] = zeroVal;
    }
}



template class Grid<float>;
template class Grid<std::vector<FluidParticle*, std::allocator<FluidParticle*> > >;
