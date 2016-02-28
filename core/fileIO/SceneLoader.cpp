//
// Created by austin on 2/26/16.
//

#include <fstream>
#include <core/geometry/Box.h>
#include "SceneLoader.h"

const char* SceneLoader::defaultScene =
        "{\n"
        "\t\"containerDim\" : {\n"
        "\t\t\"scaleX\" : 5.0,\n"
        "\t\t\"scaleY\" : 5.0,\n"
        "\t\t\"scaleZ\" : 5.0\n"
        "\t},\n"
        "\t\"particleDim\" : {\n"
        "\t\t\"boundX\" : 3.7,\n"
        "\t\t\"boundY\" : 3.7,\n"
        "\t\t\"boundZ\" : 3.7\n"
        "\t},\n"
        "\t\"particleSeparation\" : 0.1\n"
        "}";

Fluid::Solver* SceneLoader::LoadScene(const char *filepath) {
    Json::Reader reader;

    std::ifstream fileStream(filepath, std::ifstream::binary);

    Json::Value root;
    if (filepath != nullptr) {
        if (!reader.parse(fileStream, root, false)) {
            fprintf(stderr, "Failed to load json file %s!", filepath);
            exit(EXIT_FAILURE);
        }
    } else {
        return LoadScene(std::string(defaultScene));
    }

    return parseJson(root);
}

Fluid::Solver* SceneLoader::LoadScene(const std::string &jsonstring) {
    Json::Reader reader;
    Json::Value root;

    if (!reader.parse(jsonstring, root, false)) {
        fprintf(stderr, "Failed to load json string!\n %s", jsonstring);
        exit(EXIT_FAILURE);
    }

    return parseJson(root);
}

Fluid::Solver* SceneLoader::parseJson(const Json::Value &root) {
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

    Fluid::Solver* solver = new Fluid::Solver(particleSep);
    solver->setContainer(container);
    solver->addFluid(fluidObject);
    return solver;
}
