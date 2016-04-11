//
// Created by austin on 3/21/16.
//

#include "ParticlesWriter.h"
//#include <openvdb_points/tools/PointDataGrid.h>
//#include <openvdb_points/tools/PointConversion.h>
//#include <openvdb_points/tools/PointCount.h>

//using namespace openvdb::tools;

ParticlesWriter::ParticlesWriter() {
    openvdb::initialize();
//    openvdb::points::initialize();
}

ParticlesWriter::~ParticlesWriter() {

}

void ParticlesWriter::writeData(const FluidSolver* const solver, const std::string &filename) {
//    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
//    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
//    openvdb::Coord xyz(0, 0, 0);

//    for (const FluidParticle &particle : solver->_particles) {
//        xyz.reset(particle.cell.x, particle.cell.y, particle.cell.z);
//        accessor.setValue(xyz, 1.f);
//    }
//
//    openvdb::io::File file(filename);
//    openvdb::GridPtrVec grids;
//    grids.push_back(grid);
//
//    file.write(grids);
//    file.close();

//    const float voxelSize = 10.0f;
//    openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(voxelSize);
//
//    std::vector<openvdb::Vec3f> positions;
//    for (const FluidParticle &particle : solver->_particles) {
//        positions.push_back(openvdb::Vec3f(particle.pos.x, particle.pos.y, particle.pos.z));
//    }
//
//    PointDataGrid::Ptr grid = createPointDataGrid<PointDataGrid>(positions, TypedAttributeArray<openvdb::Vec3f>::attributeType(), *transform);
//
//    openvdb::io::File file("filename");
//    openvdb::GridPtrVec grids;
//    grids.push_back(grid);
//
//    file.write(grids);
//    file.close();
}
