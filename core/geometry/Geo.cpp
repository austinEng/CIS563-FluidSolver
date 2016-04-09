//
// Created by austin on 2/27/16.
//


#include "Geo.h"

bool Geo::collidesPt(const glm::vec3 &pt, float tolerance) const {
    glm::vec3 norm;
    return collidesPt(pt, norm, tolerance);
}

bool Geo::collides(const glm::vec3 &prev, const glm::vec3 &next) const {
    glm::vec3 norm;
    return collides(prev, next, norm);
}

bool Geo::collidesRay(const glm::vec3 &pt, const glm::vec3 &dir, float step) const {
    assert(fequal(glm::length(dir), 1.f));
    glm::vec3 nextPt = pt + dir * step;
    return collides(pt, nextPt);
}

bool Geo::collidesRay(const glm::vec3 &pt, const glm::vec3 &dir, glm::vec3 &normal, float step) const {
    assert(fequal(glm::length(dir), 1.f));
    glm::vec3 nextPt = pt + dir * step;
    return collides(pt, nextPt, normal);
}