#include <core/display/Window.h>
#include <core/fileIO/SceneLoader.h>

int main(int argc, char* argv[]) {
    Window* window = new Window("Fluid Solver");
    SceneLoader::LoadScene(argv[1]);

    window->initloop();
    delete window;

    return 0;
}