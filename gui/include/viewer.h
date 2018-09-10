#pragma once

#include <nanogui/screen.h>
#include <set>
#include <memory>

#include <camera.h>
#include <gl_utils.h>
#include <string>


class SceneObject; // forward declaration due to circular dependency
class ProjectGui;

typedef bool (*callback_function)(void*);

// The main class for a GUI. There can only be one of these in a given program.
class Viewer {

    public:

        // Create a new viewer object
        Viewer(std::string name="Awesome code");
        ~Viewer();

        // Start main loop for the viewer
        void startMainLoop();
        void setMainLoopCallback(callback_function func, void* payload);


        void addSceneObject(SceneObject* obj);
        void removeSceneObject(SceneObject* obj);

        // The camera object containing data about the user's view of the scene
        Camera camera;

        // A global reference to the One Viewer
        static Viewer* masterViewer;

        // The primary screen object that holds widgets for nanogui
        // Note: This object never gets deleted, but I couldn't find anyway to delete it on exit which
        // doesn't create errors with glfw in some way or other. The relevant nanogui example does the
        // same thing.
        static nanogui::Screen *mainScreen;

        // display-related member
        static GLFWwindow* mainWindow;
        static long maxLoopsPerSecond; // if lazy-drawing how many times should we loop per second (default: 100)
    
        // Global draw options
        static bool lazyDrawing; // if true, only redraw something something has changed (default: true)
        static bool showShadows;
        static bool showEdges;
        static bool flatShaded;
        static bool showDual;
        static bool pickingEnabled;
        static float wireframeOpacity;
        
        // Draw!
        void draw();

        // Dump a screenshot to disk, WITHOUT gui widgets
        void screenshot(const std::string& filename = "" );
        
        // Collection of GUIs for individual projects
        std::vector<ProjectGui*> projectGuis;

    private:

        friend class Gui;

        // Various properties
        std::string windowName;

        // Static members:
        // Note: Because the windowing infrastructure uses payload-less callbacks,
        // a few of the core data members need to be static. This limits us to a single
        // window fow now, but that should be okay. To overcome this limitation, we could
        // instead use a static list of screens/windows/etc.

        // display-related member
        std::set<SceneObject*> sceneObjects;

        // Internal callback functions
        static void cursorPositionCallback(GLFWwindow *w, double x, double y);
        static int prevMouseX;
        static int prevMouseY;
        static void mouseButtonPressCallback(GLFWwindow *w, int button, int action, int modifiers);
        static void keyPressCallback(GLFWwindow *w, int key, int scancode, int action, int mods);
        static void characterCallback(GLFWwindow *w, unsigned int codepoint);
        static void scrollCallback(GLFWwindow *w, double x, double y);
        static void framebufferSizeCallback(GLFWwindow *w, int width, int height);
        static void windowPosCallback(GLFWwindow *w, int xPos, int yPos);
        static bool anyCallbackFired;
        callback_function applicationLoopCallback = nullptr;
        void* applicationLoopCallbackPayload = nullptr;


        // Setup and initialization
        void setup();

        // Draw!
        void drawScene();

        // Call any methods that should be continuously updated
        bool animate();
    
        // Processes picked element
        void processPickedElement(const int pickedId, int x, int y);

        // dim screen by given factor
        GLProgram *dimShader;
        void dim( double amount );

};
