//
// Created by austin on 3/20/16.
//

#ifndef FLUIDSOLVER_MACGRID_H
#define FLUIDSOLVER_MACGRID_H

#include "Grid.h"
#include <core/util/flags.h>

enum CellType {
    EMPTY,
    FLUID,
    SOLID
};

template <typename T> class MACGrid : public Grid<T> {
    friend class FluidSolver;
public:

    MACGrid();
    MACGrid(const glm::vec3 &origin, const glm::vec3 &dim, float size);
    virtual ~MACGrid();

    Grid<float> _gU;
    Grid<float> _gV;
    Grid<float> _gW;
    Grid<float> _gU_old;
    Grid<float> _gV_old;
    Grid<float> _gW_old;
    Grid<float> _gP;
    Grid<int> _gType;
//    Grid<float> _gDiv;

private:
    std::vector<T> _contents;
};


#endif //FLUIDSOLVER_MACGRID_H
