#include <core/display/Window.h>
#include <core/fileIO/SceneLoader.h>
#include <core/display/painters/ParticlesPainter.h>
#include <core/display/painters/BoxPainter.h>

int main(int argc, char* argv[]) {
    Window* window = new Window("Fluid Solver");

    FluidSolver* solver = SceneLoader::LoadScene(argv[1]);

    ParticlesPainter particlesPainter(solver);
    BoxPainter boxPainter((Box *) solver->_container);

    window->addPainter(&particlesPainter);
    window->addPainter(&boxPainter);

    int framerate = 24;

    double start = glfwGetTime();
    window->initloop([&]() {
        double now = glfwGetTime();
        float duration = (float) (now - start);

        // limit solver update to 24fps
        if (duration >= 1.f / framerate) {
            start = now;
            solver->update(duration);
        }
    });
    delete window;

    return 0;
}