#pragma once

#include <nanogui/screen.h>
#include <nanogui/window.h>
#include <nanogui/layout.h>
#include <nanogui/button.h>

#include <gui.h>

#include <string>

// Abstract base class for individual project GUIs.
// Each project is responsible for implementing (i) a "show"
// method that sets up the widgets and (ii) a "revert" method
// that updates state when the mesh is reset.

class ProjectGui {

    public:
        ProjectGui(Gui* mainGui, std::string name);
        ~ProjectGui();

        void toggleVisible(void);

        virtual void revert(void) = 0;
        virtual void updatePickedElement(const int pickedId, int x, int y) {};
        virtual bool animate(void) {return false;} // return true if redraw is needed

        std::string name;

        bool getVisibility(void) const;

    protected:
        Gui* mainGui;
        nanogui::Window* window;
        bool isVisible;

        // show() is called whenever a project GUI is loaded
        // from the tool chest; it should create a window for
        // the project and initialize any necessary data (note
        // that show() may be called multiple times per session,
        // so make sure not to initialize the same data twice!)
        virtual bool show(void) = 0;

        // hide() is called whenever the window for a project
        // GUI is closed; it does not need to deallocate the
        // data for the project (which can be persistent even
        // when the window is closed), but it should remove
        // any scene objects from the main viewer scene (i.e.,
        // nothing from that project should be drawn anymore)
        virtual void hide(void);
};
