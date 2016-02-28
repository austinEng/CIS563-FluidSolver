//
// Created by austin on 2/26/16.
//

#ifndef FLUIDSOLVER_SCENELOADER_H
#define FLUIDSOLVER_SCENELOADER_H

#include <json/json.h>
#include <core/solver/FluidSolver.h>

class SceneLoader {
public:
    static Fluid::Solver* LoadScene(const char* filepath);
    static Fluid::Solver* LoadScene(const std::string &jsonstring);

    static const char * defaultScene;

private:
    static Fluid::Solver* parseJson(const Json::Value &root);
};


#endif //FLUIDSOLVER_SCENELOADER_H
