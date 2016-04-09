//
// Created by austin on 2/28/16.
//

#ifndef FLUIDSOLVER_GEOOBJECT_H
#define FLUIDSOLVER_GEOOBJECT_H

#include "Geo.h"
#include "Bound.h"

class GeoObject : public Geo {
public:
    GeoObject() { }
    virtual ~GeoObject() {}
    virtual void computeBound() = 0;
    const Bound& bound() const { return _bound; }

protected:
    Bound _bound;
};


#endif //FLUIDSOLVER_GEOOBJECT_H
