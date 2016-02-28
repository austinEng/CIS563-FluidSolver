//
// Created by austin on 2/27/16.
//

#include "Bound.h"

Bound::Bound() {}

Bound::Bound(const glm::vec3 &center, const glm::vec3 &dim) : Bound(center, dim.x, dim.y, dim.z) {}

Bound::Bound(const glm::vec3 &center, float sX, float sY, float sZ) : Bound(center.x - sX / 2.f, center.y - sY / 2.f, center.z - sZ / 2.f,
                                                                      center.x + sX/2.f, center.y + sY/2.f, center.z + sZ/2.f) {}

Bound::Bound(float cX, float cY, float cZ, const glm::vec3 &dim) : Bound(cX - dim.x, cY - dim.y, cZ - dim.z,
                                                                   cX + dim.x, cY + dim.y, cZ + dim.z) {}

Bound::Bound(float minX, float minY, float minZ, float maxX, float maxY, float maxZ)
        : _minX(minX), _minY(minY), _minZ(minZ), _maxX(maxX), _maxY(maxY), _maxZ(maxZ) {}

Bound::~Bound() {}

float Bound::minX() const { return _minX; }

float Bound::minY() const { return _minY; }

float Bound::minZ() const { return _minZ; }

float Bound::maxX() const { return _maxX; }

float Bound::maxY() const { return _maxY; }

float Bound::maxZ() const { return _maxZ; }

float Bound::width() const { return _maxX - _minX; }

float Bound::height() const { return _maxY - _minY; }

float Bound::depth() const { return _maxZ - _minZ; }

glm::vec3 Bound::dim() const { return glm::vec3(_maxX - _minX, _maxY - _minY, _maxZ - _minZ); }

glm::vec3 Bound::center() const { return glm::vec3((_minX + _maxX) / 2.f, (_minY + _maxY) / 2.f, (_minZ + _maxZ) / 2.f); }

bool Bound::contains(const glm::vec3 &pt) {
    return (pt.x >= _minX && pt.x < _maxX &&
            pt.y >= _minY && pt.y < _maxY &&
            pt.z >= _minZ && pt.z < _maxZ);
}

bool Bound::collidesPt(const glm::vec3 &pt, glm::vec3 &normal, float tolerance) {
    if (fequal(pt.x, _minX, tolerance)) {
        if (pt.y >= _minY && pt.y < _maxY && pt.z >= _minZ && pt.z < _maxZ) {
            if (pt.x < _minX) {
                normal = glm::vec3(-1.f, 0.f, 0.f);
            } else {
                normal = glm::vec3(1.f, 0.f, 0.f);
            }
            return true;
        }
    } else if (fequal(pt.x, _maxX, tolerance)) {
        if (pt.y >= _minY && pt.y < _maxY && pt.z >= _minZ && pt.z < _maxZ) {
            if (pt.x > _maxX) {
                normal = glm::vec3(1.f, 0.f, 0.f);
            } else {
                normal = glm::vec3(-1.f, 0.f, 0.f);
            }
            return true;
        }
    }
    if (fequal(pt.y, _minY, tolerance)) {
        if (pt.x >= _minX && pt.x < _maxX && pt.z >= _minZ && pt.z < _maxZ) {
            if (pt.y < _minY) {
                normal = glm::vec3(0.f, -1.f, 0.f);
            } else {
                normal = glm::vec3(0.f, 1.f, 0.f);
            }
            return true;
        }
    } else if (fequal(pt.y, _maxY, tolerance)) {
        if (pt.x >= _minX && pt.x < _maxX && pt.z >= _minZ && pt.z < _maxZ) {
            if (pt.y > _maxY) {
                normal = glm::vec3(0.f, 1.f, 0.f);
            } else {
                normal = glm::vec3(0.f, -1.f, 0.f);
            }
            return true;
        }
    }
    if (fequal(pt.z, _minZ, tolerance)) {
        if (pt.x >= _minX && pt.x < _maxX && pt.y >= _minY && pt.y < _maxY) {
            if (pt.z < _minZ) {
                normal = glm::vec3(0.f, 0.f, -1.f);
            } else {
                normal = glm::vec3(0.f, 0.f, 1.f);
            }
            return true;
        }
    } else if (fequal(pt.z, _maxZ, tolerance)) {
        if (pt.x >= _minX && pt.x < _maxX && pt.y >= _minY && pt.y < _maxY) {
            if (pt.z > _maxZ) {
                normal = glm::vec3(0.f, 0.f, 1.f);
            } else {
                normal = glm::vec3(0.f, 0.f, -1.f);
            }
            return true;
        }
    }
    return false;
}

bool Bound::collides(const glm::vec3 &prev, const glm::vec3 &next, glm::vec3 &normal) {
    if (prev.y >= _minY && prev.y < _maxY && prev.z >= _minZ && prev.z < _maxZ) {
        if (prev.x > _minX && next.x <= _minX) {        // cross minX plane
            normal = glm::vec3(1.f, 0.f, 0.f);
            return true;
        } else if (prev.x < _minX && next.x >= _minX) {
            normal = glm::vec3(-1.f, 0.f, 0.f);
            return true;
        } else if (prev.x < _maxX && next.x >= _maxX) { // cross maxX plane
            normal = glm::vec3(-1.f, 0.f, 0.f);
            return true;
        } else if (prev.x > _maxX && next.x <= _maxX) {
            normal = glm::vec3(1.f, 0.f, 0.f);
            return true;
        }
    }
    if (prev.x >= _minX && prev.x < _maxX && prev.z >= _minZ && prev.z < _maxZ) {
        if (prev.y > _minY && next.y <= _minY) {        // cross minY plane
            normal = glm::vec3(0.f, 1.f, 0.f);
            return true;
        } else if (prev.y < _minY && next.y >= _minY) {
            normal = glm::vec3(0.f, -1.f, 0.f);
            return true;
        } else if (prev.y < _maxY && next.y >= _maxY) { // cross maxY plane
            normal = glm::vec3(0.f, -1.f, 0.f);
            return true;
        } else if (prev.y > _maxY && next.y <= _maxY) {
            normal = glm::vec3(0.f, 1.f, 0.f);
            return true;
        }
    }
    if (prev.x >= _minX && prev.x < _maxX && prev.y >= _minY && prev.y < _maxY) {
        if (prev.z > _minZ && next.z <= _minZ) {        // cross minZ plane
            normal = glm::vec3(0.f, 1.f, 0.f);
            return true;
        } else if (prev.z < _minZ && next.z >= _minZ) {
            normal = glm::vec3(0.f, -1.f, 0.f);
            return true;
        } else if (prev.z < _maxZ && next.z >= _maxZ) { // cross maxZ plane
            normal = glm::vec3(0.f, -1.f, 0.f);
            return true;
        } else if (prev.z > _maxZ && next.z <= _maxZ) {
            normal = glm::vec3(0.f, 1.f, 0.f);
            return true;
        }
    }
    return false;
}