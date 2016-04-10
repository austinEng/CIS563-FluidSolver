#include <core/display/Window.h>
#include <core/fileIO/SceneLoader.h>
#include <core/fileIO/ParticlesWriter.h>
#include <core/display/painters/ParticlesPainter.h>
#include <core/display/painters/GridVectorAttributePainter.h>
#include <core/display/painters/GridScalarAttributePainter.h>
#include <core/display/painters/BoxPainter.h>
#include <ctime>

int main(int argc, char* argv[]) {
    Window* window = new Window("Fluid Solver");

    FluidSolver* solver = SceneLoader::LoadScene(argv[1]);

    ParticlesWriter particlesWriter;
    solver->init();

    particlesWriter.writeData(solver, "particles_0.vdb");

    ParticlesPainter particlesPainter(solver);
    BoxPainter boxPainter((Box *) solver->_container);
    GridVectorAttributePainter uPainter (&solver->_MAC._gU, 3.f, glm::vec3(1,0,0), glm::vec3(0.2,0,0));
    GridVectorAttributePainter vPainter (&solver->_MAC._gV, 3.f, glm::vec3(0,1,0), glm::vec3(0,0.2,0));
    GridVectorAttributePainter wPainter (&solver->_MAC._gW, 3.f, glm::vec3(0,0,1), glm::vec3(0,0,0.2));
    GridScalarAttributePainter tPainter (
            &solver->_MAC._gType, 0.f, 1.f, 0.f, 7.f, glm::vec3(0,0,0), glm::vec3(0,0,1));
    GridScalarAttributePainter posPressurePainter (
            &solver->_MAC._gP, 0.f, 100.f, 0.f, 5.f, glm::vec3(1,1,0), glm::vec3(1,0,0));
    GridScalarAttributePainter negPressurePainter (
            &solver->_MAC._gP, 0.f, -100.f, 0.f, 5.f, glm::vec3(0,1,1), glm::vec3(0,0,1));

    window->addPainter(&particlesPainter);
    window->addPainter(&boxPainter);
    window->addPainter(&uPainter);
    window->addPainter(&vPainter);
    window->addPainter(&wPainter);
//    window->addPainter(&tPainter);
    window->addPainter(&posPressurePainter);
    window->addPainter(&negPressurePainter);

    window->loadSceneCB = [](void*) {
        std::cout << "what" << std::endl;
    };
//    window->initializeTweakBar();

    int framerate = 24;
    double start = glfwGetTime();
    int frame = 0;

    float totalComputeTime = 0;

    window->initloop([&]() {
        double now = glfwGetTime();
        float duration = (float) (now - start);

        // limit solver update to 24fps
        if (duration >= 1.f / framerate) {
            start = now;
            //solver->update(duration);
            solver->update(1.f / framerate);
            totalComputeTime += glfwGetTime() - start;

            std::string filename = "particles_";
            filename.append(std::to_string(++frame));
            filename.append(".vdb");
            particlesWriter.writeData(solver, filename);
        }
    });
//    solver->update(0.1f);

    std::cout << "Average compute time: " << totalComputeTime / frame << " seconds" << std::endl;

    delete window;
    delete solver;

    return 0;
}