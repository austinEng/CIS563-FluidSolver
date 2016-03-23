//
// Created by austin on 2/28/16.
//

#ifndef FLUIDSOLVER_FLUIDSOLVER_H
#define FLUIDSOLVER_FLUIDSOLVER_H

#include <core/util/math.h>
#include <core/geometry/GeoObject.h>
#include <vector>
#include "grid/MACGrid.h"
#include "FluidParticle.h"

class FluidSolver {
    friend class ParticlesPainter;
    friend class ParticlesWriter;
public:
    FluidSolver(float particleSep, float gridSize);
    ~FluidSolver();

    void setContainer(GeoObject* container);
    void addFluid(GeoObject* fluid);
    void init();

    void projectVelocitiesToGrid();
    void transferVelocitiesToParticles();
    void gravitySolve(float step);
    void updateParticlePositions(float step);
    void resolveCollisions();
    void updateCells();

    void update(float step = 0.04166f);

    GeoObject* _container;
    MACGrid<std::vector<FluidParticle*> > _MAC;

private:
    std::vector<FluidParticle> _particles;
    float particle_radius;
    float _cell_size;
    int frame;

    template <class T> void particleAttributeToGrid(std::size_t offset, Grid<T> &grid, float radius, T zeroVal);
    template <class T> T interpolateAttribute(const glm::vec3 &pos, Grid<T> &grid);

    void iterParticles(const std::function<void(FluidParticle &particle)> &cb, bool parallel=true);

    static float g;
};


#endif //FLUIDSOLVER_FLUIDSOLVER_H
