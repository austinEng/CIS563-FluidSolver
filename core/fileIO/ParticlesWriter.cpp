//
// Created by austin on 3/21/16.
//

#include "ParticlesWriter.h"

ParticlesWriter::ParticlesWriter() {
    openvdb::initialize();
    grid = openvdb::FloatGrid::create();
}

ParticlesWriter::~ParticlesWriter() {

}

void ParticlesWriter::writeData(const FluidSolver* const solver, const std::string &filename) {
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    openvdb::Coord xyz(0, 0, 0);

    for (const FluidParticle &particle : solver->_particles) {
        xyz.reset(particle.cell.x, particle.cell.y, particle.cell.z);
        accessor.setValue(xyz, 1.f);
    }

    openvdb::io::File file(filename);
    openvdb::GridPtrVec grids;
    grids.push_back(grid);

    file.write(grids);
    file.close();
}
