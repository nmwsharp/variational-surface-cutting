#pragma once

#include <cstdlib>
#include <vector>
#include <memory>

#include <scene_object.h>
#include <gl_utils.h>
#include <viewer.h>
#include <vector3.h>
#include <shaders.h>

class SceneLines: public SceneObject {

    public:

        SceneLines(Viewer &parent);

        virtual void draw() override;
        //virtual void initialize() override;

        void setLines(std::vector<std::array<Vector3, 2>> const &lines);
    
        bool enabled = true;
        double lineRadius = 0.02;
        Vector3 lineColor = RGB_WHITE;

    private:

        // Drawing options
        GLProgram glProgram;
        unsigned int nPoints;

        // openGL things
        std::vector<Vector3> posData;
        std::vector<Vector3> normData;
        std::vector<Vector3> colorData;
        std::vector<uint3> indData;


};
