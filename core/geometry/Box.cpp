//
// Created by austin on 2/28/16.
//

#include "Box.h"

Box::Box(const glm::vec3 &center, const glm::vec3 &dim) : Box(center, dim.x, dim.y, dim.z) {}

Box::Box(const glm::vec3 &center, float sX, float sY, float sZ) : Box(center.x - sX / 2.f, center.y - sY / 2.f, center.z - sZ / 2.f,
                                                                        center.x + sX/2.f, center.y + sY/2.f, center.z + sZ/2.f) {}

Box::Box(float cX, float cY, float cZ, const glm::vec3 &dim) : Box(cX - dim.x, cY - dim.y, cZ - dim.z,
                                                                     cX + dim.x, cY + dim.y, cZ + dim.z) {}

Box::Box(float minX, float minY, float minZ, float maxX, float maxY, float maxZ) :
    Bound(minX, minY, minZ, maxX, maxY, maxZ), GeoObject() {
    computeBound();
}
