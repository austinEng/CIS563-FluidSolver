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

    window->initloop();
    delete window;

    return 0;
}