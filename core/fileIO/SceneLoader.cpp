//
// Created by austin on 2/26/16.
//

#include <fstream>
#include <core/geometry/Box.h>
#include "SceneLoader.h"
#include <core/scenes/default.h>

FluidSolver* SceneLoader::LoadScene(const char *filepath) {
    Json::Reader reader;

    std::ifstream fileStream(filepath, std::ifstream::binary);

    Json::Value root;
    if (filepath != nullptr) {
        if (!reader.parse(fileStream, root, false)) {
            fprintf(stderr, "Failed to load json file %s!", filepath);
            exit(EXIT_FAILURE);
        }
    } else {
        return LoadScene(std::string(default_scene));
    }

    return parseJson(root);
}

FluidSolver* SceneLoader::LoadScene(const std::string &jsonstring) {
    Json::Reader reader;
    Json::Value root;

    if (!reader.parse(jsonstring, root, false)) {
        fprintf(stderr, "Failed to load json string!\n %s", jsonstring);
        exit(EXIT_FAILURE);
    }

    return parseJson(root);
}

FluidSolver* SceneLoader::parseJson(const Json::Value &root) {
    Json::Value containerDim = root["containerDim"];
    Json::Value particleDim = root["particleDim"];
    Json::Value resolution = root["resolution"];

    glm::vec3 containerSize(containerDim["scale"][0].asFloat(),
                            containerDim["scale"][1].asFloat(),
                            containerDim["scale"][2].asFloat());
    glm::vec3 containerPos(containerDim["position"][0].asFloat(),
                           containerDim["position"][1].asFloat(),
                           containerDim["position"][2].asFloat());

    glm::vec3 fluidSize(particleDim["scale"][0].asFloat(),
                        particleDim["scale"][1].asFloat(),
                        particleDim["scale"][2].asFloat());
    glm::vec3 fluidPos(particleDim["position"][0].asFloat(),
                       particleDim["position"][1].asFloat(),
                       particleDim["position"][2].asFloat());

    float cellSize = std::max(std::max(containerSize.x, containerSize.y), containerSize.z) / resolution.asFloat();

    Box* container = new Box(containerPos, containerSize);
//    Box* container = new Box(containerPos, containerSize + 2.f*glm::vec3(cellSize, cellSize, cellSize));
    Box fluidObject = Box(fluidPos, fluidSize);

    FluidSolver* solver = new FluidSolver(cellSize/2, cellSize);
    solver->setContainer(container);
    solver->addFluid(fluidObject);
    return solver;
}
