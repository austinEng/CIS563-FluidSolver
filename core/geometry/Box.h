//
// Created by austin on 2/28/16.
//

#ifndef FLUIDSOLVER_BOX_H
#define FLUIDSOLVER_BOX_H

#include "GeoObject.h"

class Box : public Bound, public GeoObject {
public:
    Box(const glm::vec3 &center, const glm::vec3 &dim);
    Box(const glm::vec3 &center, float sX, float sY, float sZ);
    Box(float cX, float cY, float cZ, const glm::vec3 &dim);
    Box(float minX, float minY, float minZ, float maxX, float maxY, float maxZ);

    virtual bool contains(const glm::vec3 &pt) const {
        return _bound.contains(pt);
    }
    virtual bool collidesPt(const glm::vec3 &pt, glm::vec3 &normal, float tolerance = 0.001f) const {
        return _bound.collidesPt(pt, normal, tolerance);
    }
    virtual bool collides(const glm::vec3 &prev, const glm::vec3 &next, glm::vec3 &normal) const {
        return _bound.collides(prev, next, normal);
    }

    virtual void computeBound() {
        _bound = *this;
    }
};


#endif //FLUIDSOLVER_BOX_H
