//
// Created by austin on 3/20/16.
//

#include "Grid.h"
#include <core/solver/FluidParticle.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range3d.h>
#include <core/util/hacks.h>
#include <iostream>

template <typename T> Grid<T>::Grid() {

}

template <typename T> Grid<T>::~Grid() {

}

template <typename T> Grid<T>::Grid(const glm::vec3 &origin, const glm::vec3 &offset, const glm::vec3 &dim, float size) :
        _origin(origin),
        _offset(offset),
        _dim(dim),
        _cellSize(size),
        _countX((size_t) (std::ceil((_dim.x - _offset.x) / _cellSize))),
        _countY((size_t) (std::ceil((_dim.y - _offset.y) / _cellSize))),
        _countZ((size_t) (std::ceil((_dim.z - _offset.z) / _cellSize))) {
    _contents = std::vector<T>((unsigned long) (_countX * _countY * _countZ));
}

template <typename T> template <typename C> Grid<T>::Grid(const Grid<C> &rhs) :
        _origin(rhs._origin), 
        _offset(rhs._offset), 
        _dim(rhs._dim), 
        _cellSize(rhs._cellSize),
        _countX(rhs._countX),
        _countY(rhs._countY),
        _countZ(rhs._countZ) {
    _contents = std::vector<T>((unsigned long) (_countX * _countY * _countZ));
}

template <typename T> T& Grid<T>::operator()(std::size_t idx) {
    return _contents[idx];
}

template <typename T> const T& Grid<T>::operator()(std::size_t idx) const {
    return _contents[idx];
}

template <typename T> T& Grid<T>::operator()(std::size_t i, std::size_t j, std::size_t k) {
    return _contents[k*_countX*_countY + j*_countX + i];
}

template <typename T> const T& Grid<T>::operator()(std::size_t i, std::size_t j, std::size_t k) const {
    return _contents[k*_countX*_countY + j*_countX + i];
}


template <typename T> T& Grid<T>::atIdx(std::size_t i, std::size_t j, std::size_t k) {
    return _contents[k*_countX*_countY + j*_countX + i];
}

template <typename T> const T& Grid<T>::atIdx(std::size_t i, std::size_t j, std::size_t k) const {
    return _contents[k*_countX*_countY + j*_countX + i];
}


template <typename T> T& Grid<T>::operator()(const glm::ivec3 &idx) {
    return _contents[idx.z*_countX*_countY + idx.y*_countX + idx.x];
}

template <typename T> const T& Grid<T>::operator()(const glm::ivec3 &idx) const {
    return _contents[idx.z*_countX*_countY + idx.y*_countX + idx.x];;
}

template <typename T> T& Grid<T>::atIdx(const glm::ivec3 &idx) {
    return _contents[idx.z*_countX*_countY + idx.y*_countX + idx.x];
}

template <typename T> const T& Grid<T>::atIdx(const glm::ivec3 &idx) const {
    return _contents[idx.z*_countX*_countY + idx.y*_countX + idx.x];
}


template <typename T> T& Grid<T>::at(float x, float y, float z) {
    return at(glm::vec3(x, y, z));
}

template <typename T> const T& Grid<T>::at(float x, float y, float z) const {
    return at(glm::vec3(x, y, z));;
}

template <typename T> T& Grid<T>::at(const glm::vec3 &pos) {
    glm::ivec3 indices = indexOf(pos);
    return this->operator()((size_t) indices.x, (size_t) indices.y, (size_t) indices.z);
}

template <typename T> const T& Grid<T>::at(const glm::vec3 &pos) const {
    glm::ivec3 indices = indexOf(pos);
    return this->operator()((size_t) indices.x, (size_t) indices.y, (size_t) indices.z);
}

template <typename T> glm::ivec3 Grid<T>::indexOf(const glm::vec3 &pos) const {
    glm::vec3 indices = (pos - _offset - _origin) / _cellSize;
    int i = (int) indices.x;
    int j = (int) indices.y;
    int k = (int) indices.z;
    //if (i >= _countX ) i = -1;
    //if (j >= _countY ) j = -1;
    //if (k >= _countZ ) k = -1;
    return glm::clamp(glm::ivec3(i, j, k), glm::ivec3(0,0,0), glm::ivec3(_countX-1, _countY-1, _countZ-1));
}

template <typename T> void Grid<T>::indexOf(const glm::vec3 &pos, size_t &i, size_t &j, size_t &k) const {
    glm::vec3 indices = (pos - _offset - _origin) / _cellSize;
    i = (size_t) ((indices.x <= _countX) * indices.x + (indices.x > _countX) * _countX); // clamp at countX
    j = (size_t) ((indices.y <= _countY) * indices.y + (indices.y > _countY) * _countY); // clamp at countY
    k = (size_t) ((indices.x <= _countZ) * indices.z + (indices.z > _countZ) * _countZ); // clamp at countZ
    i = (i > 0) * i;
    j = (j > 0) * j;
    k = (k > 0) * k;
}


template <typename T> glm::vec3 Grid<T>::positionOf(const glm::ivec3 &idx) const {
    return glm::vec3(idx.x * _cellSize, idx.y * _cellSize, idx.z * _cellSize) + _offset + _origin;
}

template <typename T> glm::vec3 Grid<T>::positionOf(size_t i, size_t j, size_t k) const {
    return glm::vec3(i * _cellSize, j * _cellSize, k * _cellSize) + _offset + _origin;
}

template <typename T> glm::vec3 Grid<T>::fractionalIndexOf(const glm::vec3 &pos) const {
    return glm::clamp((pos - _offset - _origin) / _cellSize, glm::vec3(0,0,0), glm::vec3(_countX, _countY, _countZ));
}

template <typename T> glm::ivec3 Grid<T>::toIJK(const std::size_t index) const {
    int i = (int) (index % _countZ);
    int j = (int) ((index / _countZ) % _countY);
    int k = (int) (index / (_countY * _countZ));
    return glm::ivec3(i,j,k);
}

template <typename T> std::size_t Grid<T>::fromIJK(const std::size_t i, const std::size_t j, const std::size_t k) const {
    return (size_t) (k * _countX * _countY + j * _countX + i);
}

template <typename T> std::size_t Grid<T>::fromIJK(const glm::ivec3 &ijk) const {
    return (size_t) (ijk.z * _countX * _countY + ijk.y * _countX + ijk.x);
}

template <typename T> void Grid<T>::iterate(const std::function<void(size_t i, size_t j, size_t k)> &cb, bool parallel) {
#ifdef USETBB
    if (parallel) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, _contents.size()), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                glm::ivec3 ijk = toIJK(i);
                cb(ijk.x, ijk.y, ijk.z);
            }
        });
    } else {
        for (size_t idx = 0; idx < _contents.size(); idx++) {
            glm::ivec3 ijk = toIJK(idx);
            cb(ijk.x, ijk.y, ijk.z);
        }
    }
#else
    for (size_t idx = 0; idx < _contents.size(); idx++) {
        glm::ivec3 ijk = toIJK(idx);
        cb(ijk.x, ijk.y, ijk.z);
    }
#endif
}

template <typename T> void Grid<T>::iterateNeighborhood(size_t i, size_t j, size_t k, size_t r, const std::function<void(size_t i, size_t j, size_t k)> &cb, bool parallel) {
    size_t si = MATHIFELSE(i - r, 0, i == 0);
    size_t sj = MATHIFELSE(j - r, 0, j == 0);
    size_t sk = MATHIFELSE(k - r, 0, k == 0);
    size_t ei = MATHIFELSE(i + r, _countX-1, i + r >= _countX);
    size_t ej = MATHIFELSE(j + r, _countY-1, j + r >= _countY);
    size_t ek = MATHIFELSE(k + r, _countZ-1, k + r >= _countZ);

#ifdef USETBB
    if (parallel) {
        tbb::parallel_for(tbb::blocked_range3d<size_t>(si,ei,sj,ej,sk,ek), [&](const tbb::blocked_range3d<size_t> &r) {
            for(size_t i=r.pages().begin(), i_end=r.pages().end(); i<i_end; i++){
                for (size_t j=r.rows().begin(), j_end=r.rows().end(); j<j_end; j++){
                    for (size_t k=r.cols().begin(), k_end=r.cols().end(); k<k_end; k++){
                        std::cout << i << "," << j << "," << k << std::endl;
                        cb(i,j,k);
                    }
                }
            } 
        });
    } else {
        for (size_t i = si; i <= ei; i++) {
            for (size_t j = sj; j <= ej; j++) {
                for (size_t k = sk; k <= ek; k++) {
                    cb(i,j,k);
                }
            }
        }
    }
#else
    for (size_t i = si; i <= ei; i++) {
        for (size_t j = sj; j <= ej; j++) {
            for (size_t k = sk; k <= ek; k++) {
                cb(i,j,k);
            }
        }
    }
#endif
}

template <typename T> void Grid<T>::getNeighboorhood(size_t i, size_t j, size_t k, size_t r, size_t &si, size_t &ei, size_t &sj, size_t &ej, size_t &sk, size_t &ek) {
    si = MATHIFELSE(i - r, 0, i == 0);
    sj = MATHIFELSE(j - r, 0, j == 0);
    sk = MATHIFELSE(k - r, 0, k == 0);
    ei = std::min(i+r+1, _countX); //MATHIFELSE(i + r, _countX, i + r >= _countX);
    ej = std::min(j+r+1, _countY); //MATHIFELSE(j + r, _countY, j + r >= _countY);
    ek = std::min(k+r+1, _countZ); //MATHIFELSE(k + r, _countZ, k + r >= _countZ);
}


template <typename T> void Grid<T>::clear(const T &zeroVal) {
#ifdef USETBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, _contents.size()), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            _contents[i] = zeroVal;
        }
    });
#else
    for (size_t i = 0; i < _contents.size(); i++) {
        _contents[i] = zeroVal;
    }
#endif
}

template <typename T> bool Grid<T>::checkIdx(size_t i, size_t j, size_t k) const {
    return i >= 0 && i < _countX &&
            j >= 0 && j < _countY &&
            k >= 0 && k < _countZ;
}
template <typename T> bool Grid<T>::checkIdx(const glm::ivec3 &idx) const {
    return checkIdx((size_t) idx.x, (size_t) idx.y, (size_t) idx.z);
}


template class Grid<float>;
template class Grid<std::vector<FluidParticle*, std::allocator<FluidParticle*> > >;