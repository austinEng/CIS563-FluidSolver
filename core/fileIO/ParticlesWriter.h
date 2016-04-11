//
// Created by austin on 3/21/16.
//

#ifndef FLUIDSOLVER_PARTICLESWRITER_H
#define FLUIDSOLVER_PARTICLESWRITER_H

#include <openvdb/openvdb.h>
//#include <openvdb_points/openvdb.h>
#include <core/solver/FluidSolver.h>

class ParticlesWriter {

public:
    ParticlesWriter();
    ~ParticlesWriter();

    void writeData(const FluidSolver* const solver, const std::string &filename);
};


#endif //FLUIDSOLVER_PARTICLESWRITER_H
