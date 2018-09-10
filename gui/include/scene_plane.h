#pragma once

#include <cstdlib>
#include <vector>
#include <memory>

// Nanogui
#include <nanogui/button.h>
#include <nanogui/label.h>
#include <nanogui/window.h>
#include <nanogui/textbox.h>
#include <nanogui/layout.h>
#include <nanogui/graph.h>
#include <nanogui/checkbox.h>
#include <nanogui/slider.h>

#include <scene_object.h>
#include <gl_utils.h>
#include <colormaps.h>
#include <viewer.h>
#include <vector3.h>
#include <shaders.h>

class ScenePlane: public SceneObject {

    public:
        ScenePlane(Viewer &parent);

        virtual void draw() override;
        virtual void drawLight() override;

        // Drawing options
        bool showPlane;

        // Specify shift and two basis vectors; scale facilitates shading
        void setGeometry( double offset, Vector3 e1, Vector3 e2, Vector3 center, double scale = 1. );

    private:
        double scale;
        GLProgram shader;
        float intensity;

        void drawPlane();

};
