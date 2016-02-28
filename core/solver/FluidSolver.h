//
// Created by austin on 2/28/16.
//

#ifndef FLUIDSOLVER_FLUIDSOLVER_H
#define FLUIDSOLVER_FLUIDSOLVER_H

#include <core/util/math.h>
#include <core/geometry/GeoObject.h>
#include <vector>

namespace Fluid {
    struct Particle {
        glm::vec3 pos;
        glm::vec3 vel;
        glm::vec3 col;

        Particle() {
            pos = glm::vec3(0);
            vel = glm::vec3(0);
            col = glm::vec3(0.5f, 0.5f, 1.f);
        }
    };

    class Solver {
    public:
        Solver(float particleSep);
        ~Solver();

        void setContainer(GeoObject* container);
        void addFluid(GeoObject* fluid);

        void update(float step = 0.04166f);

    private:
        GeoObject* _container;
        std::vector<Particle*> _particles;
        float particle_radius;

        void presolve(float step);
        void solve(float step);
        void postsolve(float step);

        static float g;
    };
}


#endif //FLUIDSOLVER_FLUIDSOLVER_H
