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
    Json::Value particleSeparation = root["particleSeparation"];

    float containerX = containerDim["scaleX"].asFloat();
    float containerY = containerDim["scaleY"].asFloat();
    float containerZ = containerDim["scaleZ"].asFloat();

    float particleX = particleDim["boundX"].asFloat();
    float particleY = particleDim["boundY"].asFloat();
    float particleZ = particleDim["boundZ"].asFloat();

    float particleSep = particleSeparation.asFloat();

    glm::vec3 origin;
    Box* container = new Box(origin, containerX, containerY, containerZ);
    Box* fluidObject = new Box(origin, particleX, particleY, particleZ);

    FluidSolver* solver = new FluidSolver(particleSep);
    solver->setContainer(container);
    solver->addFluid(fluidObject);
    return solver;
}
