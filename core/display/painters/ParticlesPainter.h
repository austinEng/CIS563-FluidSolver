//
// Created by austin on 2/29/16.
//

#ifndef FLUIDSOLVER_PARTICLESPAINTER_H
#define FLUIDSOLVER_PARTICLESPAINTER_H

#include "Painter.h"
#include <core/solver/FluidSolver.h>

class ParticlesPainter : public Painter {
public:
    ParticlesPainter(FluidSolver* solver, float ptSize = 3.f);
    virtual void draw() const;

private:
    unsigned int MAX_PARTICLES = 10000;
    GLuint particle_buffer;

    GLint attrPos;
    GLint attrVel;
    GLint attrCol;

    GLfloat _ptSize;
    std::vector<FluidParticle>* _particles;
};


#endif //FLUIDSOLVER_PARTICLESPAINTER_H
