//
// Created by austin on 2/27/16.
//

#ifndef FLUIDSOLVER_GEO_H
#define FLUIDSOLVER_GEO_H

#include <core/util/math.h>

class Geo {
public:
    Geo() { }
    virtual ~Geo() { }

    virtual bool contains(const glm::vec3 &pt) = 0;

    virtual bool collidesPt(const glm::vec3 &pt, float tolerance = 0.001f);

    virtual bool collidesPt(const glm::vec3 &pt, glm::vec3 &normal, float tolerance = 0.001f) = 0;

    virtual bool collides(const glm::vec3 &prev, const glm::vec3 &next);

    virtual bool collides(const glm::vec3 &prev, const glm::vec3 &next, glm::vec3 &normal) = 0;

    virtual bool collidesRay(const glm::vec3 &pt, const glm::vec3 &dir, float step = 0.001f);

    virtual bool collidesRay(const glm::vec3 &pt, const glm::vec3 &dir, glm::vec3 &normal, float step = 0.001f);
};


#endif //FLUIDSOLVER_GEO_H
