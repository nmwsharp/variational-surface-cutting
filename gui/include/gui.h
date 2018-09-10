#pragma once

#include <viewer.h>
#include <scene_mesh.h>
#include <geometry.h>
#include <scene_plane.h>
#include <scene_points.h>
#include <string>
#include <set>

// Class which encapsulates the data for the GUI
//
// Note that this is somewhat redundant with the Viewer class; both encapsulate GUI information
// and there should only be one of them in the program. They are separated to keep things reasonably
// sized: the Viewer handles low level rendering and windowing; the Gui constructs the actual user
// interface filled with wonderful buttons and text boxes that is the primary interface for this codebase.
// One good reason to keep these classes separate is that Viewer contains broadly useful tools for setting
// up an openGL/nanogui window, whereas Gui is strictly for creating this particular interface.
class Gui {

    public:
        
        // Constructor
        Gui();

        // Main data
        Viewer* viewer = nullptr;
        
        HalfedgeMesh* mesh = nullptr;
        Geometry<Euclidean>* geometry = nullptr;
        
        SceneMesh* sceneMesh = nullptr; // scene mesh currently shown in viewer
        SceneMesh* defaultSceneMesh = nullptr; // fallback in case no application wants to display a custom scene mesh

        ScenePlane* scenePlane = nullptr; // ground plane

        // The primary panel holding display options
        nanogui::Window *primaryWindow;
        nanogui::Label* meshVerticesLabel;
        nanogui::CheckBox* meshDualCheckbox;
        nanogui::CheckBox* meshEdgeCheckbox;
        nanogui::Slider* meshEdgeOpacitySlider;
        nanogui::Label* meshEdgeOpacityLabel;

        // Utility function for easy error reporting
        // (Currently shows a pop-up window, but this message
        // could also be routed to a console, file, ...)
        void reportError(const std::string& errorMessage);

        // Things relating to submodules
        nanogui::Window *toolChestWindow;

        // Functions
        void loadMesh(std::string meshFilename);
        void exportMesh(std::string meshFilename);
        void replaceMainSceneMesh(SceneMesh* newSceneMesh);
        void replaceMainMeshAndDeleteOld(SceneMesh* newSceneMesh);
        void revertMesh();
        void addGroundPlane();
        void updateGroundPlane();
        void loadPrimaryGui();
        void loadProjectGuis();
        void buildToolChest();

    protected:

        // Path to most recently loaded mesh (for revert)
        std::string currentMeshPath;
};
