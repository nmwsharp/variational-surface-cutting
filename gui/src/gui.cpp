#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>

// Nanogui #include <nanogui/screen.h>
#include <nanogui/button.h>
#include <nanogui/label.h>
#include <nanogui/window.h>
#include <nanogui/textbox.h>
#include <nanogui/layout.h>
#include <nanogui/graph.h>
#include <nanogui/slider.h>

// Includes from this project (and core)
#include <gui.h>
#include <halfedge_mesh.h>
#include <project_gui.h>
#include <timing.h>
#include <polygon_soup_mesh.h>
#include <utilities.h>
#include <meshio.h>

#include <cuts_gui.h>

Gui::Gui() : currentMeshPath("") {
    
    // Create a viewer
    viewer = new Viewer("Grand Central Gui");
    viewer->camera.dist = 5.0;


    // === GUI stuff ===
    loadPrimaryGui();
    loadProjectGuis();
    viewer->mainScreen->performLayout(); // lay out all of the newly created elements

    // Show a ground plane for the scene
    addGroundPlane();

}

// Load a mesh in to the program
void Gui::loadMesh(std::string meshFilename) {

    std::cout << "Loading mesh " << meshFilename << std::endl;
    
    // Delete the old mesh, if there is one
    if(mesh != nullptr) {
        viewer->removeSceneObject(sceneMesh);
        delete sceneMesh;
        sceneMesh = nullptr;
        delete mesh;
        mesh = nullptr;
    }

    // Parse in the mesh
    PolygonSoupMesh inputMesh(meshFilename);
    mesh = new HalfedgeMesh(inputMesh, geometry);
    
    // Add the new mesh to the viewer
    SceneMesh* newSceneMesh = new SceneMesh(*viewer, geometry);
    newSceneMesh->resetCameraForData();
    defaultSceneMesh = newSceneMesh;

    // Keep track of the path to the mesh, in cae the user wants to revert
    currentMeshPath = meshFilename;
    
    meshVerticesLabel->setCaption("Mesh: " + std::to_string(mesh->nVertices()) + " vertices");
    
    // Adjust the ground plane height according to model size
    updateGroundPlane();

    replaceMainSceneMesh(newSceneMesh);
}
        
void Gui::replaceMainSceneMesh(SceneMesh* newSceneMesh) {

    if( newSceneMesh == nullptr ) {
        newSceneMesh = defaultSceneMesh;
        defaultSceneMesh->update();
    }
    
    viewer->removeSceneObject(sceneMesh);
    viewer->addSceneObject(newSceneMesh);
    sceneMesh = newSceneMesh;
}


void Gui::replaceMainMeshAndDeleteOld(SceneMesh* newSceneMesh) {

    // Clear out the old viz mesh
    SceneMesh* origMesh = sceneMesh;
    if(origMesh != nullptr) {
        viewer->removeSceneObject(origMesh);
        delete origMesh;
    }

    sceneMesh = newSceneMesh;
    defaultSceneMesh = newSceneMesh;
    geometry = sceneMesh->getGeometry();
    mesh = geometry->getMesh();
    
    // Set the nVerts label
    meshVerticesLabel->setCaption("Mesh: " + std::to_string(mesh->nVertices()) + " vertices");
    
    // Adjust the ground plane height according to model size
    updateGroundPlane();
}

// Export a mesh to disk
void Gui::exportMesh(std::string meshFilename) {

    std::cout << "Exporting mesh " << meshFilename << std::endl;
    
    // Make sure there's a mesh to export!
    if(mesh == nullptr) {
       reportError( "No mesh loaded." );
       return;
    }

    WavefrontOBJ::write( meshFilename, *geometry );
}

void Gui::revertMesh()
{
   if( currentMeshPath != "" ) {
      loadMesh(currentMeshPath);

      for(auto project : viewer->projectGuis)
      {
         project->revert();
      }
   }
}

void Gui::addGroundPlane()
{
    // Add the new mesh to the viewer
    scenePlane = new ScenePlane(*viewer);

    if( mesh != nullptr )
    {
       Vector3 c = geometry->center();
       Vector3 e = geometry->extent();
       double offset = c.y - e.y;
       scenePlane->setGeometry( offset, Vector3{0.,0.,1.}, Vector3{1.,0.,0.}, c, e.y );
    }
    else
    {
       scenePlane->setGeometry( -1.5, Vector3{0.,0.,1.}, Vector3{1.,0.,0.}, Vector3 {0,0,0} );
    }
}

void Gui::updateGroundPlane()
{
   if(scenePlane) {
      delete scenePlane;
   };
   addGroundPlane();
}

void Gui::loadPrimaryGui() {

    // === Create the primary gui panel
    primaryWindow = new nanogui::Window((viewer->mainScreen), "Options");
    primaryWindow->setPosition(Eigen::Vector2i(15, 15));
    primaryWindow->setLayout(new nanogui::GroupLayout());


    // === Add a button to load a mesh
    nanogui::Button *bLoadMesh = new nanogui::Button(primaryWindow, "Load mesh");
    bLoadMesh->setCallback([this] {
        std::cout << "Selecting mesh to load...  ";
        
        std::string meshFilename = nanogui::file_dialog({{"obj", "Wavefront OBJ file"}}, false);
        if(meshFilename == "") {
            std::cout << "...no file selected" << std::endl;
            return;
        }
        std::cout << "...selected mesh file " << meshFilename << std::endl;
       
        // Try to load the mesh. If something goes wrong, catch the exception instead of exiting the app.
        try {
            loadMesh(meshFilename);
        } catch (...) {
            std::cout << "Mesh loading failed." << std::endl;
        }
    });

    // === Add a button to export the current mesh mesh
    nanogui::Button *bExportMesh = new nanogui::Button(primaryWindow, "Export mesh");
    bExportMesh->setCallback([this] {
        std::cout << "Selecting mesh to export...  ";
        
        std::string meshFilename = nanogui::file_dialog({{"obj", "Wavefront OBJ file"}}, true);
        if(meshFilename == "") {
            std::cout << "...no file selected" << std::endl;
            return;
        }
        std::cout << "...selected mesh file " << meshFilename << std::endl;
       
        // Try to load the mesh. If something goes wrong, catch the exception instead of exiting the app.
        try {
            exportMesh(meshFilename);
        } catch (...) {
            std::cout << "Mesh export failed." << std::endl;
        }
    });

    // === Add a button to revert the current mesh
    nanogui::Button *bRevertMesh = new nanogui::Button(primaryWindow, "Revert mesh");
    bRevertMesh->setCallback([this] {
       revertMesh();
    });
    
    // === Add a button to take a screenshot (sans GUI elements)
    nanogui::Button *bScreenshot = new nanogui::Button(primaryWindow, "Take screenshot");
    bScreenshot->setCallback([this] {
       viewer->screenshot();
    });


    long nVerts = (mesh == nullptr) ? 0 : mesh->nVertices();
    meshVerticesLabel = new nanogui::Label(primaryWindow, "Mesh:  " + std::to_string(nVerts) + " vertices", "sans-bold");
    nanogui::Widget *meshWidget = new nanogui::Widget(primaryWindow);
    nanogui::GridLayout *meshLayout = new nanogui::GridLayout(nanogui::Orientation::Horizontal, 2, nanogui::Alignment::Middle);
    meshLayout->setColAlignment({ nanogui::Alignment::Maximum, nanogui::Alignment::Fill });
    meshLayout->setSpacing(0, 10);
    meshWidget->setLayout(meshLayout);


//   	meshDualCheckbox = new nanogui::CheckBox(primaryWindow, "Draw dual",
// 		[this](bool state) { 
// 			if(state) {
// 				std::cout << "GUI: Showing mesh dual " << std::endl; 
// 			} else {
// 				std::cout << "GUI: Hiding mesh dual " << std::endl; 
// 			}
//                         viewer->showDual = state;
// 		}
// 	);
// 	meshDualCheckbox->setChecked(false);


    // === Add options for the ground plane
  	nanogui::CheckBox *cbPlane = new nanogui::CheckBox(primaryWindow, "Draw ground plane",
		[this](bool state) { 
			if(state) {
				std::cout << "GUI: Showing ground plane " << std::endl; 
			} else {
				std::cout << "GUI: Hiding ground plane " << std::endl; 
			}
            if(scenePlane != nullptr) {
                scenePlane->showPlane = state;
            }
		}
	);
	cbPlane->setChecked(true);
        
    // === Add options for shadows
  	nanogui::CheckBox *cbShadow = new nanogui::CheckBox(primaryWindow, "Draw shadows",
		[this](bool state) { 
			if(state) {
				std::cout << "GUI: Showing shadows " << std::endl; 
			} else {
				std::cout << "GUI: Hiding shadows " << std::endl; 
			}
                viewer->showShadows = state;
		}
	);
	cbShadow->setChecked(viewer->showShadows);

    // === Add options for mesh edges
        meshEdgeOpacityLabel = new nanogui::Label(primaryWindow,"Wireframe");

  	meshEdgeCheckbox = new nanogui::CheckBox(primaryWindow, "Show",
		[this](bool state) { 
                   viewer->showEdges = state;
                   viewer->flatShaded = viewer->showEdges; // These parameters don't have to be coupled, but we'll keep it this way for now
                   if( sceneMesh != nullptr ) {
                      sceneMesh->update();
                   }
		}
	);
	meshEdgeCheckbox->setChecked(viewer->showEdges);

        meshEdgeOpacitySlider = new nanogui::Slider(primaryWindow);
        meshEdgeOpacitySlider->setValue(viewer->wireframeOpacity);
        meshEdgeOpacitySlider->setCallback(
              [this](float value) {
                 viewer->wireframeOpacity = value;
                 }
              );
}


void Gui::buildToolChest() {
   
   toolChestWindow = new nanogui::Window((viewer->mainScreen), "Tool Chest");
   toolChestWindow->setPosition(Eigen::Vector2i(200,15));
   toolChestWindow->setLayout(new nanogui::GroupLayout());

   for(auto project : viewer->projectGuis) {
      nanogui::Button *bToggleProject = new nanogui::Button(toolChestWindow, project->name);
      bToggleProject->setFlags( nanogui::Button::ToggleButton );
      bToggleProject->setChangeCallback([project,bToggleProject](bool){
            project->toggleVisible();
            bToggleProject->setPushed(project->getVisibility());
            });
   }
}

// Perform GUI intialization for each external gui
void Gui::loadProjectGuis() {

    viewer->projectGuis.push_back( new CutsGui(this,"Variational Cuts") );

    buildToolChest();
}
