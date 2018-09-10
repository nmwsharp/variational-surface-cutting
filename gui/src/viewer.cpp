#include<viewer.h>

#include <gui.h>
#include <image.h>

#include "project_gui.h"

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <thread>

#include <scene_object.h>
#include <gl_utils.h>
#include <GLFW/glfw3.h>
#include <colors.h>
#include "timing.h"

// Nanogui
#include <nanogui/screen.h>
#include <nanogui/button.h>
#include <nanogui/label.h>
#include <nanogui/window.h>
#include <nanogui/layout.h>

#include <shaders/dim_shader.h>
#define MAX_PICKED_ID 16777215

// Initialize statics
nanogui::Screen* Viewer::mainScreen = new nanogui::Screen(); // Note: can't find a way to safely delete this. See note in header.
int Viewer::prevMouseX = 0;
int Viewer::prevMouseY = 0;
long Viewer::maxLoopsPerSecond = 100;
bool Viewer::lazyDrawing = true;
bool Viewer::anyCallbackFired = false;
bool Viewer::showShadows = false;
bool Viewer::showEdges = false;
bool Viewer::flatShaded = false;
bool Viewer::showDual = false;
bool Viewer::pickingEnabled = false;
float Viewer::wireframeOpacity = .5;
GLFWwindow* Viewer::mainWindow;
Viewer* Viewer::masterViewer;

Viewer::Viewer(std::string name)
    : windowName(name), dimShader( NULL )
{
    // Peform setup tasks
    setup();
}

Viewer::~Viewer()
{
   if( dimShader ) {
      delete dimShader;
   }
}

void Viewer::addSceneObject(SceneObject* obj) {
    sceneObjects.insert(obj);
}
void Viewer::removeSceneObject(SceneObject* obj) {
    for(auto o : sceneObjects) {
        if(o == obj) {
            sceneObjects.erase(o);
            return;
        }
    }
}

//  === Callbacks ===

// Small callback function for GLFW errors
void error_print_callback(int error, const char* description) {
    std::cerr << "GLFW emitted error: " << description << std::endl;
}


void Viewer::cursorPositionCallback(GLFWwindow *w, double x, double y) {
    anyCallbackFired = true;

    bool nanoGuiHit = mainScreen->cursorPosCallbackEvent(x,y);

    // Do other stuff
    if(!nanoGuiHit) {
        // Pass to camera if appropriate
        if(glfwGetMouseButton(mainWindow, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
            int xSize, ySize;
            glfwGetWindowSize(mainWindow, &xSize, &ySize);

            double oldX = (double) prevMouseX / xSize;
            double newX = (double) x / xSize;
            double oldY = (double) prevMouseY / ySize;
            double newY = (double) y / ySize;

            // Process the movement as a translation (rather than a rotation) if either of
            // the shift keys are currently held
            int leftShiftState = glfwGetKey(mainWindow, GLFW_KEY_LEFT_SHIFT);
            int rightShiftState = glfwGetKey(mainWindow, GLFW_KEY_RIGHT_SHIFT);
            bool isRotating = !(leftShiftState == GLFW_PRESS || rightShiftState == GLFW_PRESS);

            //Viewer::masterViewer->camera.mouseDragEvent(delX, delY, isRotating);
            Viewer::masterViewer->camera.mouseDragEvent(oldX, oldY, newX, newY, isRotating);
        }
    }

    prevMouseX = x;
    prevMouseY = y;
}

void Viewer::mouseButtonPressCallback(GLFWwindow *w, int button, int action, int modifiers) {
    anyCallbackFired = true;

    // Propagate the event to the nanogui UI
    bool nanoGuiHit = mainScreen->mouseButtonCallbackEvent(button, action, modifiers);
    if(nanoGuiHit) {
        return;
    }

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) pickingEnabled = true;
}

void Viewer::keyPressCallback(GLFWwindow *w, int key, int scancode, int action, int mods) {
    anyCallbackFired = true;

    // Propagate the event to the nanogui UI
    bool nanoGuiHit = mainScreen->keyCallbackEvent(key, scancode, action, mods);
    if (nanoGuiHit) {
        return;
    }

    // Application level commands (just quit for now)
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(mainWindow, true);
    }
}
       
// Note: The distinction between this and the keyPressCallback is that that this handles
// textual input (like typing words), while the keyPressCallback handles modifiers and special
// keys. 
void Viewer::characterCallback(GLFWwindow *w, unsigned int codepoint) {
    anyCallbackFired = true;

    // Propagate the event to the nanogui UI
    bool nanoGuiHit = mainScreen->keyboardCharacterEvent(codepoint);
    if(nanoGuiHit) {
        return;
    }

    // Do other stuff
}

void Viewer::scrollCallback(GLFWwindow *w, double x, double y) {
    anyCallbackFired = true;

    // Propagate the event to the nanogui UI
    bool nanoGuiHit = mainScreen->scrollCallbackEvent(x, y);
    if(nanoGuiHit) {
        return;
    }

    // On some setups, shift flips the scroll direction, so take the max scrolling in any direction
    double maxScroll = x;
    if(std::abs(y) > std::abs(x)) {
        maxScroll = y;
    }

    // Pass camera commands to the camera
    if(maxScroll != 0.0) {

        int leftShiftState = glfwGetKey(mainWindow, GLFW_KEY_LEFT_SHIFT);
        int rightShiftState = glfwGetKey(mainWindow, GLFW_KEY_RIGHT_SHIFT);
        bool scrollClipPlane = (leftShiftState == GLFW_PRESS || rightShiftState == GLFW_PRESS);
    
        if(scrollClipPlane) {
            Viewer::masterViewer->camera.mouseScrollEvent(maxScroll, scrollClipPlane); // shift could be flipping scroll direction
        } else {
            Viewer::masterViewer->camera.mouseScrollEvent(y, scrollClipPlane);
        }
    }
}

void Viewer::framebufferSizeCallback(GLFWwindow *w, int width, int height) {
    anyCallbackFired = true;

    // Pass camera commands to the camera
    Viewer::masterViewer->camera.setWindowSize(width, height);

    // Propagate the event to the nanogui UI
    mainScreen->resizeCallbackEvent(width, height);

    masterViewer->draw();
}

void Viewer::windowPosCallback(GLFWwindow *w, int xPos, int yPos) {
    anyCallbackFired = true;
    masterViewer->draw();
}


void Viewer::setup() {

    Viewer::masterViewer = this;

    std::cout << "Creating window and initalizing openGL..." << std::endl;
    
    // Attempt to initialize GLFW
    glfwSetErrorCallback(error_print_callback);
    if (!glfwInit()) {
        std::cerr << "GLFW failed to initialize...";
        exit(EXIT_FAILURE);
    }

    // Set hints for the GLFW window
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwWindowHint(GLFW_SAMPLES, 0);
    glfwWindowHint(GLFW_RED_BITS, 8);
    glfwWindowHint(GLFW_GREEN_BITS, 8);
    glfwWindowHint(GLFW_BLUE_BITS, 8);
    glfwWindowHint(GLFW_ALPHA_BITS, 8);
    glfwWindowHint(GLFW_STENCIL_BITS, 8);
    glfwWindowHint(GLFW_DEPTH_BITS, 24);
    glfwWindowHint(GLFW_VISIBLE, GL_FALSE);
    //glfwWindowHint(GLFW_RESIZABLE, GL_FALSE); 
   

    // Create a main GLFW window
    const int initWindowWidth = 1500;
    const int initWindowHeight = 900;
    mainWindow = glfwCreateWindow(initWindowWidth, initWindowHeight, windowName.c_str(), nullptr, nullptr);
    if (!mainWindow) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    // Activate the context
    glfwMakeContextCurrent(mainWindow);
    glfwSwapInterval(1);
    

    // Register callbacks
    glfwSetCursorPosCallback(mainWindow, Viewer::cursorPositionCallback);
    glfwSetMouseButtonCallback(mainWindow, Viewer::mouseButtonPressCallback);
    glfwSetKeyCallback(mainWindow, Viewer::keyPressCallback);
    glfwSetCharCallback(mainWindow, Viewer::characterCallback);
    glfwSetScrollCallback(mainWindow, Viewer::scrollCallback);
    glfwSetFramebufferSizeCallback(mainWindow, Viewer::framebufferSizeCallback);
    glfwSetWindowPosCallback(mainWindow, Viewer::windowPosCallback);

    
    // Load openGL functions (using GLAD)
#ifndef __APPLE__
    if(!gladLoadGL()) {
        std::cerr << "ERROR: Failed to load openGL using GLAD" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Loaded openGL version: " <<  glGetString(GL_VERSION) << std::endl;
#endif


#ifdef __APPLE__
    // Hack to classify the process as interactive
    glfwPollEvents();
#endif


    // Initialize the screen
    glfwMakeContextCurrent(mainWindow);
    Viewer::mainScreen->initialize(mainWindow, false);

    // Make sure the camera knows about the window size
    camera.setWindowSize(initWindowWidth, initWindowHeight);


    // initialize shaders (must be done after creation of GL context)
    GLProgram::initCommonShaders(); // common utility functions available to all shaders
    dimShader = new GLProgram(&DIM_VERT_SHADER, &DIM_GEOM_SHADER, &DIM_FRAG_SHADER, DrawMode::Points);

    // Allocate a dummy point for drawing a full-screen quad in the dimmer shader.
    std::vector<Vector3> positions( 1 );
    positions[0] = Vector3{ 0., 0., 0. };
    dimShader->setAttribute( "a_position", positions );

    std::cout << "  ...window creation and openGL setup complete." << std::endl << std::endl;
}


// Specify an application-defined command which gets called on every iteration
void Viewer::setMainLoopCallback(callback_function func, void* payload) {
    applicationLoopCallback = func;
    applicationLoopCallbackPayload = payload;
}

// Begin the main display loop for the program. Kicks off a new thread (after completing
// initializiation tasks) and returns while that thread runs.
void Viewer::startMainLoop() {
    
    // Make sure the window is visible
    glfwShowWindow(mainWindow);


    // Layout the nanogui widgets which may have been added between initialization
    // and starting the main loop
    mainScreen->performLayout();
    mainScreen->setVisible(true);


    // Run the window in a loop
    bool redrawNeeded = true;
    while (!glfwWindowShouldClose(mainWindow))
    {
        auto loopStartTime = NOW;

        // Poll for input events, flagging for a redraw if anyhing happened
        anyCallbackFired = false;
        glfwPollEvents();
        redrawNeeded |= anyCallbackFired;


        // Call the application callback, if one is registered
        if(applicationLoopCallback != nullptr) {
            bool callbackRedrawNeeded = applicationLoopCallback(applicationLoopCallbackPayload);
            redrawNeeded |= callbackRedrawNeeded;
        }

        // Call animate loops for all projects 
        bool animateRedrawNeeded = animate();
        redrawNeeded |= animateRedrawNeeded;

        // Draw, if needed
        if(redrawNeeded) {
            draw();
        }

        // Make sure we don't loop more than X times per second
        // to avoid busy waiting while not drawing
        if(lazyDrawing) {
            auto currTime = NOW;
            long microsecPerLoop = 1000000/maxLoopsPerSecond;
            while(std::chrono::duration_cast<std::chrono::microseconds>(currTime - loopStartTime).count() < microsecPerLoop) {
                std::chrono::milliseconds timespan(1);
                std::this_thread::sleep_for(timespan);
                currTime = NOW;
            }
    
            redrawNeeded = false;
        }
    }


    // Shutting down glfw here leads to errors as stack-allocated nanogui objects get destructed later.
    // Until we need something from the glfw shutdown, let's just not do it.
    //glfwDestroyWindow(mainWindow);
    //glfwTerminate();
}

bool Viewer::animate() {

   bool anyRedraw = false;
   for( auto g : projectGuis )
   {
      anyRedraw |= g->animate();
   }

   return anyRedraw;
}

void Viewer::processPickedElement(const int pickedId, int x, int y)
{
    for (auto gui : masterViewer->projectGuis) {
        gui->updatePickedElement(pickedId, x, y);
    }
}

void Viewer::drawScene() {
   
    // Get window parametesr
    glfwMakeContextCurrent(mainWindow);
    int width, height;
    int xSize, ySize;
    glfwGetWindowSize(mainWindow, &xSize, &ySize);
    glfwGetFramebufferSize(mainWindow, &width, &height);

    // Clear the current buffer
    // NOTE: Need to be careful with window size vs. framebuffer size for hidpi rendering
    glViewport(0, 0, width, height);
    
    // Pick
    if (pickingEnabled &&
        prevMouseX >= 0 && prevMouseX <= xSize &&
        prevMouseY >= 0 && prevMouseY <= ySize) {
        glClearColor(1.0, 1.0, 1.0, 1.0);
        glClearDepth(1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        for (auto obj: sceneObjects) {
            obj->drawPick();
        }
        
        glFlush();
        glFinish();
        
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        
        unsigned char data[4];
        // glReadPixels(prevMouseX*2, height - prevMouseY*2, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data);
        glReadPixels(prevMouseX, height - prevMouseY, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data);
        
        // Convert color to ID
        int pickedId = data[0] + data[1]*256 + data[2]*256*256;
        if (pickedId != MAX_PICKED_ID) processPickedElement(pickedId, prevMouseX, prevMouseY);
        
        pickingEnabled = false;
    }
    
    Vector3 bg = getBackgroundColor();
    glClearColor( bg.x, bg.y, bg.z, 1. );
    glClearDepth( 1. );
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    
    if( showShadows )
    {
       glEnable(GL_DEPTH_TEST);
       glDepthFunc(GL_LEQUAL);

       // Draw objects that will be affected by shadows
       for(auto obj : sceneObjects) {
          if(obj->receiveShadow) {
             obj->draw();
          }
       }

       // Darken all shadowed objects
       glDepthFunc( GL_NOTEQUAL ); // don't darken the background
       glDepthMask(0); // dimmer shouldn't modify Z-buffer
       dim( .5 );

       // Now draw shadow volumes into stencil buffer
       glColorMask(0, 0, 0, 0);
       glDepthMask(0);
       glDepthFunc(GL_LESS);
       glEnable( GL_DEPTH_CLAMP );
       glEnable(GL_STENCIL_TEST);
       glStencilFunc(GL_ALWAYS, 0, ~0);
       glStencilOpSeparate(  GL_BACK, GL_KEEP, GL_INCR_WRAP, GL_KEEP);
       glStencilOpSeparate( GL_FRONT, GL_KEEP, GL_DECR_WRAP, GL_KEEP);
       glEnable( GL_POLYGON_OFFSET_FILL );
       glPolygonOffset( 0.1 , 5. );
       for(auto obj : sceneObjects) {
          obj->drawShadow();
       }
       glDisable( GL_POLYGON_OFFSET_FILL );
       glColorMask(1, 1, 1, 1);
       glDepthMask(1);
       glDisable( GL_DEPTH_CLAMP );

       // Draw lit geometry
       glDepthFunc(GL_LEQUAL);
       glStencilFunc(GL_EQUAL, 0, ~0);
       glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
       glEnable( GL_BLEND );
       glBlendFunc( GL_ONE, GL_ONE );
       for(auto obj : sceneObjects) {
          obj->drawLight();
       }
       glDisable(GL_BLEND);
       glDepthFunc(GL_LESS);
       glDisable(GL_STENCIL_TEST);
       
       // Finally, draw objects that were not shadow receivers
       for(auto obj : sceneObjects) {
          if(obj->receiveShadow == false) {
             obj->draw();
          }
       }
    }
    else
    {
       glEnable(GL_DEPTH_TEST);
       glDepthFunc(GL_LESS);
       for(auto obj : sceneObjects) {
          obj->draw();
       }
    }
}

void Viewer::screenshot(const std::string& filename) {

   drawScene();

   GLint viewport[4];
   glGetIntegerv( GL_VIEWPORT, viewport );
   int w = viewport[2];
   int h = viewport[3];

   Image image( w, h );
   glReadPixels( 0, 0, w, h, GL_BGR, GL_FLOAT, &image(0,0) );

   if( filename != "" )
   {
      // If a specific filename was specified, use that
      image.write( filename.c_str() );
   }
   else
   {
      // Otherwise, fall back to dumping to a numbered sequence

      static int nScreenshots = 0;

      std::stringstream ss;
      ss << "screenshot";
      ss << std::setw(10) << std::setfill('0') << nScreenshots;
      ss << ".tga";

      image.write( ss.str().c_str() );

      nScreenshots++;
   }
}

void Viewer::draw() {

   // Draw 3D geometry
   drawScene();

   // Draw the nanogui gui
   mainScreen->drawWidgets();

   glfwSwapBuffers(mainWindow);
}

void Viewer::dim(double amount) {
   glEnable( GL_BLEND );
   glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
   dimShader->setUniform( "u_dimAmount", amount );
   dimShader->draw();
   glDisable( GL_BLEND );
}


#include <nanogui/messagedialog.h>
void Gui::reportError(const std::string& errorMessage)
{
   nanogui::MessageDialog* errorDialog =
      new nanogui::MessageDialog((viewer->mainScreen),
         nanogui::MessageDialog::Type::Warning,
         "Error",
         errorMessage,
         "I'm sorry");
}
