//
// Created by austin on 2/27/16.
//

#ifndef FLUIDSOLVER_BOUND_H
#define FLUIDSOLVER_BOUND_H

#include "Geo.h"

class Bound : public Geo {
public:
    Bound();
    Bound(const glm::vec3 &center, const glm::vec3 &dim);
    Bound(const glm::vec3 &center, float sX, float sY, float sZ);
    Bound(float cX, float cY, float cZ, const glm::vec3 &dim);
    Bound(float minX, float minY, float minZ, float maxX, float maxY, float maxZ);
    virtual ~Bound();

    float minX() const;
    float minY() const;
    float minZ() const;
    float maxX() const;
    float maxY() const;
    float maxZ() const;
    float width() const;
    float height() const;
    float depth() const;
    glm::vec3 dim() const;
    glm::vec3 center() const;

    virtual bool contains(const glm::vec3 &pt);
    virtual bool collidesPt(const glm::vec3 &pt, glm::vec3 &normal, float tolerance = 0.001f);
    virtual bool collides(const glm::vec3 &prev, const glm::vec3 &next, glm::vec3 &normal);

private:
    float _minX;
    float _minY;
    float _minZ;
    float _maxX;
    float _maxY;
    float _maxZ;
};


#endif //FLUIDSOLVER_BOUND_H
