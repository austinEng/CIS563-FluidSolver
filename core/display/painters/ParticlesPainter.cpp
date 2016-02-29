//
// Created by austin on 2/29/16.
//

#include "ParticlesPainter.h"
#include <core/display/shaders/particle.frag.h>
#include <core/display/shaders/particle.vert.h>

ParticlesPainter::ParticlesPainter(FluidSolver* solver, float ptSize) : _ptSize(ptSize), _particles(&solver->_particles) {
    GLuint particleVert = compileShader(particle_vert, GL_VERTEX_SHADER);
    GLuint particleFrag = compileShader(particle_frag, GL_FRAGMENT_SHADER);

    std::vector<GLuint> programs = {particleVert, particleFrag};
    prog = makeProgram(programs);

    attrPos = glGetAttribLocation(prog, "v_pos");
    attrVel = glGetAttribLocation(prog, "v_vel");
    attrCol = glGetAttribLocation(prog, "v_col");

    glGenBuffers(1, &particle_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particle_buffer);
    glBufferData(GL_ARRAY_BUFFER, MAX_PARTICLES * sizeof(glm::vec3), NULL, GL_STREAM_DRAW);
}

void ParticlesPainter::draw() const {
    if (_particles != nullptr) {
        glUseProgram(prog);

        glBindBuffer(GL_ARRAY_BUFFER, particle_buffer);
        glBufferData(GL_ARRAY_BUFFER, MAX_PARTICLES * sizeof(FluidParticle), NULL, GL_STREAM_DRAW);
        glBufferSubData(GL_ARRAY_BUFFER, 0, MAX_PARTICLES * sizeof(FluidParticle), &((*_particles)[0]));

        glEnableVertexAttribArray(attrPos);
        glVertexAttribPointer(attrPos, 3, GL_FLOAT, GL_FALSE, sizeof(FluidParticle), (void*)offsetof(FluidParticle, pos));

        glEnableVertexAttribArray(attrVel);
        glVertexAttribPointer(attrVel, 3, GL_FLOAT, GL_FALSE, sizeof(FluidParticle), (void*)offsetof(FluidParticle, vel));

        glEnableVertexAttribArray(attrCol);
        glVertexAttribPointer(attrCol, 3, GL_FLOAT, GL_FALSE, sizeof(FluidParticle), (void*)offsetof(FluidParticle, col));

        glPointSize(_ptSize);
        glDrawArrays(GL_POINTS, 0, (GLsizei) _particles->size());

        glDisableVertexAttribArray(attrPos);
        glDisableVertexAttribArray(attrVel);
        glDisableVertexAttribArray(attrCol);
    }
}