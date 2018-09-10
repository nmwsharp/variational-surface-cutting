#pragma once

#include "project_gui.h"

#include "eulerian_shape_optimizer.h"

#include "scene_mesh_cuts.h"

class CutsGui : public ProjectGui {


    public:
        CutsGui(Gui* mainGui, std::string name);
        ~CutsGui();
        
        virtual void revert(void) override;

    protected:

        virtual bool show(void) override;
        virtual void hide(void) override;
        virtual void updatePickedElement(const int pickedId, int x, int y) override;

        // Use a second window
        nanogui::Window* paramWindow = nullptr;

        // Implicit: these ptrs are non-null iff the object exists
        SceneMeshCuts* sceneMeshCuts = nullptr;

        // Cut optimizers
        EulerianShapeOptimizer* shapeOpt = nullptr;

        // Function pointers for callbacks
        std::function<void(void)> *updateVisualization;
        std::function<void(void)> *updateParameters;
   
        HalfedgeMesh* mesh = nullptr;
        Geometry<Euclidean>* geometry = nullptr;
};
