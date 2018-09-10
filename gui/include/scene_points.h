#pragma once

#include <cstdlib>
#include <vector>
#include <memory>

#include <scene_object.h>
#include <gl_utils.h>
#include <viewer.h>
#include <vector3.h>
#include <shaders.h>

// A general base class for something that can be drawn in the scene
class ScenePoints : public SceneObject {

    public:

        ScenePoints(Viewer &parent);
        std::string name = "unnamed particle group";

        virtual void draw() override;
        //virtual void initialize() override;

        void setPositions(std::vector<Vector3> const &pos);
    
        bool enabled = true;
        int nThetaSphere = 8;
        int nPhiSphere = 5;
        float sphereRad = 0.02;
        Vector3 pointColor = RGB_WHITE;

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
