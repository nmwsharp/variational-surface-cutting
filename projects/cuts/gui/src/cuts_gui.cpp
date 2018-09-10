#include "cuts_gui.h"

#include <nanogui/button.h>
#include <nanogui/combobox.h>
#include <nanogui/layout.h>
#include <nanogui/screen.h>
#include <nanogui/textbox.h>
#include <nanogui/popupbutton.h>
#include <nanogui/popup.h>
#include <nanogui/widget.h>
#include <nanogui/window.h>

#include "meshio.h"

#include "sdf_initialization.h"
#include "detect_symmetry.h"
#include "generate_geometry.h"
#include "region_management.h"
#include "read_scalar_function_from_texcoord.h"

using std::cout; using std::cerr; using std::endl;

CutsGui::CutsGui(Gui* mainGui, std::string name)
: ProjectGui(mainGui, name)
{}

// Helpers
namespace {
    void initializeTextBox(nanogui::TextBox* t, double initVal) {
        t->setEditable(true);
        t->setFixedSize(Eigen::Vector2i(100, 20));
        std::string initValStr = std::to_string(initVal);
        t->setValue(initValStr);
        t->setDefaultValue(initValStr);
        t->setFontSize(16);
        // t->setFormat("[-]?[0-9]*\\.?[0-9]+") // apparently slightly old C++ compilers don't support the <regex> that this uses
        t->setFormat("");
    }
    
    void initializeTextBox(nanogui::TextBox* t, int initVal) {
        t->setEditable(true);
        t->setFixedSize(Eigen::Vector2i(100, 20));
        std::string initValStr = std::to_string(initVal);
        t->setValue(initValStr);
        t->setDefaultValue(initValStr);
        t->setFontSize(16);
        //t->setFormat("[-]?[0-9]*\\.?[0-9]+") // apparently slightly old C++ compilers don't support the <regex> that this uses
        t->setFormat("");
    }
    
    void tryParseDouble(nanogui::TextBox* t, double& x) {
        std::string str = t->value();
        try {
            x = std::stod(str);
        } catch (...) {
            std::cout << "Error: Could not parse string to a double. Value is: " << str << std::endl;
            t->setValue(t->defaultValue());
        }
    }
    
    void tryParseInt(nanogui::TextBox* t, int& x) {
        std::string str = t->value();
        try {
            x = std::stoi(str);
        } catch (...) {
            std::cout << "Error: Could not parse string to an int. Value is: " << str << std::endl;
            t->setValue(t->defaultValue());
        }
    }


    // helper class for term modules
    class TermWidget {

        public:
        
        nanogui::Widget* widget;
        nanogui::Label* titleLabel;
        nanogui::TextBox* weightBox;
        nanogui::CheckBox* localCB;
        std::function<void(void)> *pushParameters;

        TermWidget(nanogui::Widget* parent, std::string name, bool haveLocalCB=false) {
            
            new nanogui::Label(parent, ""); // spacer

            // Create the large, colored, centered title label
            nanogui::Widget* labelWidget = new nanogui::Widget(parent);
            nanogui::BoxLayout* labelLayout = new nanogui::BoxLayout(nanogui::Orientation::Vertical);
            labelWidget->setLayout(labelLayout);
            titleLabel = new nanogui::Label(labelWidget, name, "sans-bold");
            titleLabel->setFontSize(22);

            // Create a widget to hold all the fields for the term
            widget = new nanogui::Widget(parent);
            nanogui::GridLayout* layout = new nanogui::GridLayout(nanogui::Orientation::Horizontal, 2,
                                                nanogui::Alignment::Middle, 0, 5);
            layout->setColAlignment({ nanogui::Alignment::Maximum, nanogui::Alignment::Fill });
            layout->setSpacing(0, 10);
            widget->setLayout(layout);

            new nanogui::Label(widget, "weight:", "sans");
            
            weightBox = new nanogui::TextBox(widget);
            initializeTextBox(weightBox, 0.0);
            weightBox->setCallback([&](const std::string &newVal) {
                updateStatus();
                return true;
            });

            if(haveLocalCB) {
                new nanogui::Label(widget, "Use local scale:", "sans");
                localCB = new nanogui::CheckBox(widget, "");
            }

        };

        void updateStatus() {

            double weight = 0.0;
            tryParseDouble(weightBox, weight);

            if(weight == 0.0) {
                titleLabel->setColor(nanogui::Color(169, 50, 38, 255));
            } else {
                titleLabel->setColor(nanogui::Color(35, 155, 86, 255));
            }

        }
    };

  // Save ply with distortion as color
  void saveDistortionPLY(EulerianShapeOptimizer* shapeOpt, std::string filename) {

      cout << "Saving distortion color ply to " << filename << endl;

      Geometry<Euclidean>* cutGeom = shapeOpt->getCustomGeometry();
      HalfedgeMesh* cutMesh = cutGeom->getMesh();
      VertexData<double> distortion = shapeOpt->getDistortion();

      // == Options
      bool useScaleFactor = true;

      //Vector2 limits{1.0, 4.0};
      //Vector2 limits{-1.3, 1.3};
      Vector2 limits{-.7, .5};

      const Colormap& cm = CM_COOLWARM;

      // Print max and min, for debugging
      double minU =  std::numeric_limits<double>::infinity();
      double maxU = -std::numeric_limits<double>::infinity();
      for(VertexPtr v : cutMesh->vertices()) {
        double u = distortion[v];
        minU = std::min(minU, u);
        maxU = std::max(maxU, u);
      }
      cout << "u ranges from [" << minU << " , " << maxU << "]" << endl;

      // Apply colormap
      VertexData<Vector3> colors(cutMesh);
      for(VertexPtr v : cutMesh->vertices()) {

          // Plot u
          if(useScaleFactor)  {
            double stretch = std::abs(distortion[v]);
            double val = (stretch - limits[0]) / (limits[1] - limits[0]);
            Vector3 color = cm.getValue(val);
            colors[v] = color;
          } 
          // Plot e^2u
          else {
            double lengthFactor = std::abs(std::exp(2.0*distortion[v]));
            double stretch = std::max(lengthFactor, 1.0 / lengthFactor); // either "shrinking by a factor of X", or "growing by a factor of X", whichever is larger
            double val = (stretch - limits[0]) / (limits[1] - limits[0]);
            Vector3 color = cm.getValue(val);
            colors[v] = color;
          }
      }

      // Save
      PLY::write(filename, *cutGeom, colors); 
  }
}

bool CutsGui::show(void) {

    if(window != nullptr)
    {
        delete window;
        window = nullptr;
    }
    
    /*
    if(mainGui->mesh == nullptr) {
        mainGui->reportError("No mesh loaded.");
        return false; // failed to show
    }
    mesh = mainGui->mesh;
    geometry = mainGui->geometry;
    */
    
    // EXPERIMENT: Load symmetric icosahedron if no mesh specified
    SymmetryResult initSym;
    initSym.symmetryFound = false;
    if(mainGui->mesh == nullptr) {
        geometry = generateSymmetricIcosahedron(50000, initSym);
        // geometry = generateSymmetricIcosahedron(5000, initSym);
        mesh = geometry->getMesh();
        SceneMesh* newSceneMesh = new SceneMesh(*(mainGui->viewer), geometry);
        newSceneMesh->resetCameraForData();
        mainGui->replaceMainMeshAndDeleteOld(newSceneMesh);
    } else {
        mesh = mainGui->mesh;
        geometry = mainGui->geometry;
    }
    
    if(!mainGui->mesh->isSimplicial()) {
        mainGui->reportError("Requires pure triangle mesh.");
        return false; // failed to show
    }

    // Create the optimizer
    shapeOpt = new EulerianShapeOptimizer(geometry);
    if(initSym.symmetryFound) { // from experiment above
        shapeOpt->symmetry = initSym;
        shapeOpt->hasSymmetry = true;
    }

    // Create the cut visualizer
    if(sceneMeshCuts != nullptr) delete sceneMeshCuts;
    sceneMeshCuts = new SceneMeshCuts(*(mainGui->viewer), geometry);
    mainGui->replaceMainSceneMesh( sceneMeshCuts );
    mainGui->sceneMesh = sceneMeshCuts;
    

    // The main window pane which holds all control related things
    window = new nanogui::Window((mainGui->viewer->mainScreen), "Cuts [Control]");
    window->setPosition(Eigen::Vector2i(15, 200));
    window->setLayout(new nanogui::GroupLayout());
    
    
    // The window pane which holds all parameter related things
    paramWindow = new nanogui::Window((mainGui->viewer->mainScreen), "Cuts [Parameters]");
    paramWindow->setPosition(Eigen::Vector2i(200, 200));
    paramWindow->setLayout(new nanogui::GroupLayout());

    // Gui elements that will be filled out below
    // (need the pointers here so we can use them in labmdas)
    static nanogui::ComboBox* poissonMethodCombo;
    static nanogui::PopupButton *vizPopupPanelButton;
    static nanogui::ComboBox* vectorVisualizeCombo;
    static nanogui::CheckBox *showLinesCB, *showExtraLinesCB; 
    static nanogui::CheckBox *updateRenderCB, *useLineSearchCB, *makeMovieCB, *rotateCB, *saveEveryCB; 
    
    static nanogui::TextBox *stepSizeBox, *alphaLengthBox, *alphaOptBox, *nStepsBox, *kInspectBox;
    static nanogui::TextBox *patchZeroAreaBox;
    static nanogui::Label* nPatchLabel;
    static std::string vizName = "Region Labels";

    static std::vector<TermWidget*> termWidgets; // Note: this list leaks memory right now
    termWidgets.clear();
    
    updateParameters = new std::function<void()>([&]() {

        // Push step size
        tryParseDouble(stepSizeBox, shapeOpt->stepSizeParam);

        // Push term values
        for(TermWidget* widget : termWidgets) {
            (*widget->pushParameters)();
        }
    });
        

    updateVisualization = new std::function<void()>([&]() {

        cout << "Updating visualization..." << endl;

        // Get the region index for the moethds that need it
        int kInspect;
        tryParseInt(kInspectBox, kInspect);
        // kInspect = clamp(kInspect, 0, K_REGIONS-1);

        Geometry<Euclidean>* origGeom = mainGui->geometry;

        // Visualize appropriate surface data
        if(vizName == "Distortion") {
            Geometry<Euclidean>* cutGeom = shapeOpt->getCustomGeometry();
            sceneMeshCuts->setCurrentMesh(cutGeom);
            sceneMeshCuts->showScalar(shapeOpt->getDistortion(), true);

        } else if (vizName == "Curvature") {
            sceneMeshCuts->setCurrentMesh(origGeom);
            sceneMeshCuts->showScalar(shapeOpt->faceCurvature, true);

        } else if (vizName == "Phi") {
            sceneMeshCuts->setCurrentMesh(origGeom);
            sceneMeshCuts->showScalar(shapeOpt->getPhiK(kInspect), true);
        
        } else if (vizName == "Cut Mesh") {
            Geometry<Euclidean>* cutGeom = shapeOpt->getCustomGeometry();
            sceneMeshCuts->setCurrentMesh(cutGeom);
            sceneMeshCuts->showNothing();
        
        } else if (vizName == "Cut Mesh Param Coarse") {
            CornerData<Vector2> paramCoords;
            Geometry<Euclidean>* cutGeom = shapeOpt->getCoarseFlattenedGeometry(paramCoords); 
            sceneMeshCuts->setCurrentMesh(cutGeom);
            sceneMeshCuts->showRegionChecker(paramCoords);
        
        } else if (vizName == "Cut Mesh Param") {
            CornerData<Vector2> paramCoords;
            Geometry<Euclidean>* cutGeodesicGeom = shapeOpt->getGeodesicFlattenedGeometry(paramCoords); 
            sceneMeshCuts->setCurrentMesh(cutGeodesicGeom);
            sceneMeshCuts->showRegionChecker(paramCoords);

        } else if (vizName == "Developable Approximation") {
            Geometry<Euclidean>* devGeom = shapeOpt->getExtrinsicDevelopableGeometry(); 
            sceneMeshCuts->setCurrentMesh(devGeom);
            sceneMeshCuts->showNothing();

        } else if (vizName == "Mesh-k") {
            Geometry<Euclidean>* cutRegionGeom = shapeOpt->getPatchGeometry(kInspect);
            sceneMeshCuts->setCurrentMesh(cutRegionGeom);
            sceneMeshCuts->showNothing();
        
        } else if (vizName == "Shape Motion-k") {
            Geometry<Euclidean>* cutGeom = shapeOpt->getCustomGeometry();
            sceneMeshCuts->setCurrentMesh(cutGeom);
            sceneMeshCuts->showScalar(shapeOpt->getRegionShapeOptimizerMotion(kInspect), true);
        
        } else if (vizName == "Shape Motion") {
            sceneMeshCuts->setCurrentMesh(origGeom);
            sceneMeshCuts->showScalar(shapeOpt->getShapeOptimizerMotion());
        
        } else if (vizName == "Patch Area") {
            Geometry<Euclidean>* cutGeom = shapeOpt->getCustomGeometry();
            sceneMeshCuts->setCurrentMesh(cutGeom);
            sceneMeshCuts->showScalar(shapeOpt->getPatchArea());
        
        } else if (vizName == "Patch Normal Deviation") {
            Geometry<Euclidean>* cutGeom = shapeOpt->getCustomGeometry();
            sceneMeshCuts->setCurrentMesh(cutGeom);
            sceneMeshCuts->showScalar(shapeOpt->getPatchNormalDeviation());

        } else if (vizName == "Region Labels") {
            sceneMeshCuts->setCurrentMesh(origGeom);
            sceneMeshCuts->showRegionLabels(shapeOpt->getRegionValues()); 

        } else if (vizName == "Boundary Distance") {
            sceneMeshCuts->setCurrentMesh(origGeom);
            sceneMeshCuts->showScalar(shapeOpt->getBoundaryDistance()); 
        
        } else if (vizName == "Scaled Distortion") {
            Geometry<Euclidean>* cutGeom = shapeOpt->getCustomGeometry();
            sceneMeshCuts->setCurrentMesh(cutGeom);
            sceneMeshCuts->showScalar(shapeOpt->getScaledDistortion(), true); 

        } else if (vizName == "Feature Size") {
            sceneMeshCuts->setCurrentMesh(origGeom);
            sceneMeshCuts->showScalar(shapeOpt->featureSize); 
        
        } else if (vizName == "Visibility") {
            sceneMeshCuts->setCurrentMesh(origGeom);
            sceneMeshCuts->showScalar(shapeOpt->getVisibility()); 
        
        } else if (vizName == "Symmetry") {
            sceneMeshCuts->setCurrentMesh(origGeom);
            if(shapeOpt->hasSymmetry) {
                sceneMeshCuts->showSymmetry(shapeOpt->symmetry);
            } else {
                sceneMeshCuts->showNothing();
            }
        } else if (vizName == "DEBUG") {
            // Geometry<Euclidean>* cutGeom = shapeOpt->getCustomGeometry();
            // sceneMeshCuts->setCurrentMesh(cutGeom);
            // sceneMeshCuts->showScalar(shapeOpt->getDebug());
            
            sceneMeshCuts->setCurrentMesh(origGeom);
            sceneMeshCuts->showScalar(shapeOpt->getDebug());
        } else if (vizName == "None") {
            sceneMeshCuts->setCurrentMesh(origGeom);
            sceneMeshCuts->showNothing();

        } else {
            throw std::invalid_argument("Unrecognized mesh viz name [" + vizName + "]");
        }


        // Manage the lines toggles
        bool cutLinesToggle = showLinesCB->checked();
        if(cutLinesToggle) {
            sceneMeshCuts->setLines(shapeOpt->getBoundaryLines()); 
        }
        sceneMeshCuts->toggleLines(cutLinesToggle);
        
        bool extraCutLinesToggle = showExtraLinesCB->checked();
        if(extraCutLinesToggle) {
            sceneMeshCuts->setLinesAlternate(shapeOpt->getExtraCutLines()); 
        }
        sceneMeshCuts->toggleLinesAlternate(extraCutLinesToggle);


        // Manage widgets
        for(TermWidget* widget : termWidgets) {
            widget->updateStatus();
        }

        // Update # patches textbox
        shapeOpt->ensureHaveConnectedComponents(); // TODO for now, always update
        if(shapeOpt->nConnectedComponents > 0) {
            nPatchLabel->setCaption(std::to_string(shapeOpt->nConnectedComponents));
        } else {
            nPatchLabel->setCaption("not computed");
        }

    });
    
    new nanogui::Label(window, "Initialization", "sans-bold");
    
    // Popup for initialization
    nanogui::PopupButton *popupInitializeButton = new nanogui::PopupButton(window, "Initialize", ENTYPO_ICON_DATABASE);
    nanogui::Popup *popup = popupInitializeButton->popup();
    popup->setLayout(new nanogui::GroupLayout());

    { // Random
        nanogui::Button *bInitRandom = new nanogui::Button(popup, "Random");
        bInitRandom->setCallback([this] {
            VertexData<LabelVec> randPhi = randomMSDF(geometry);
            shapeOpt->setState(randPhi);
            shapeOpt->iIter = 0;
            (*updateVisualization)();
        });
    }
    
    { // Normal clustering 
        nanogui::Button *bInitNormal = new nanogui::Button(popup, "Normal clustering");
        bInitNormal->setCallback([this] {
            VertexData<LabelVec> normalPhi = normalClusterMSDF(geometry, PI/3.0);
            shapeOpt->setState(normalPhi);
            shapeOpt->iIter = 0;
            (*updateVisualization)();
        });
    }
    
    { // Distance clustering 
        nanogui::Button *bInitNormal = new nanogui::Button(popup, "Distance clustering");
        bInitNormal->setCallback([this] {
            VertexData<LabelVec> phi = distanceClusterMSDF(geometry, 30);
            shapeOpt->setState(phi);
            shapeOpt->iIter = 0;
            (*updateVisualization)();
        });
    }
    
    { // Mean curvature thresh 
        nanogui::Button *bInitNormal = new nanogui::Button(popup, "Mean curvature level set");
        bInitNormal->setCallback([this] {
            VertexData<LabelVec> phi = meanCurvatureLevelSetMSDF(geometry, 1.0);
            shapeOpt->setState(phi);
            shapeOpt->iIter = 0;
            (*updateVisualization)();
        });
    }
    
    { // Volleyball 
        nanogui::Button *bInitNormal = new nanogui::Button(popup, "Volleyball");
        bInitNormal->setCallback([this] {
            VertexData<LabelVec> phi = volleyballMSDF(geometry);
            shapeOpt->setState(phi);
            shapeOpt->iIter = 0;
            (*updateVisualization)();
        });
    }
    
    
    { // Single vertex 
        nanogui::Button *bInitRandom = new nanogui::Button(popup, "Single vertex");
        bInitRandom->setCallback([this] {
            // VertexData<LabelVec> phi = singleVertex(geometry);
            VertexData<LabelVec> phi = singleVertex(geometry, 32958);
            // VertexData<LabelVec> phi = singleVertex(geometry, 4474);
            // VertexData<LabelVec> phi = singleVertex(geometry, 18453); // botijo_big
            // VertexData<LabelVec> phi = singleVertex(geometry, 139617); // car hood
            shapeOpt->setState(phi);
            shapeOpt->iIter = 0;
            (*updateVisualization)();
        });
    }
    
    { // Cube
        nanogui::Button *cubeInit = new nanogui::Button(popup, "Cube");
        cubeInit->setCallback([this] {
            // VertexData<LabelVec> phi = cubeMSDF(geometry);
            VertexData<LabelVec> phi = cubeCircleMSDF(geometry);
            shapeOpt->setState(phi);
            shapeOpt->iIter = 0;
            (*updateVisualization)();
        });
    }

    { // (sorta) D-charts 
        nanogui::Button *bInitDCharts = new nanogui::Button(popup, "D-Charts");
        bInitDCharts->setCallback([this] {
            VertexData<LabelVec> dChartsPhi = dChartsMSDF(geometry);
            shapeOpt->setState(dChartsPhi);
            shapeOpt->iIter = 0;
            (*updateVisualization)();
        });
    }
    
    { // Left-right
        nanogui::Button *bInitNormal = new nanogui::Button(popup, "X-axis split");
        bInitNormal->setCallback([this] {
            VertexData<LabelVec> splitPhi = xAxisSplitMSDF(geometry);
            shapeOpt->setState(splitPhi);
            shapeOpt->iIter = 0;
            (*updateVisualization)();
        });
    }

    { // Sphere around .5,.5,.5
        nanogui::Button *bInitNormal = new nanogui::Button(popup, ".5-center sphere");
        bInitNormal->setCallback([this] {
            VertexData<LabelVec> splitPhi = sphereAtMSDF(geometry, Vector3{0.5,0.5,0.0}, 0.25);
            shapeOpt->setState(splitPhi);
            shapeOpt->iIter = 0;
            (*updateVisualization)();
        });
    }

    { // Scalar function from file
        nanogui::Button *bLoadScalar = new nanogui::Button(popup, "Load 'curvature' from file");
        bLoadScalar->setCallback([this] {
            std::cout << "Selecting curvature function load...  ";
            
            std::string filename = nanogui::file_dialog({{"obj", "Wavefront OBJ"}}, false);
            if(filename == "") {
                std::cout << "...no file selected" << std::endl;
                return;
            }
            std::cout << "...selected curvature file " << filename << std::endl;
        
            // Try to load the state. If something goes wrong, catch the exception instead of exiting.
            try {

                // Parse
                VertexData<double> f = readScalarFromTex(mesh, filename);

                // Set
                shapeOpt->setFakeFaceCuvature(f); 

            } catch (const std::exception &exc) {
                std::cout << "State loading failed." << std::endl;
                mainGui->reportError(exc.what());
            } catch (...) {
                std::cout << "State loading failed." << std::endl;
            }
            
            (*updateVisualization)();
        });
    }
    
    { // Scalar function from z 
        nanogui::Button *bLoadScalar = new nanogui::Button(popup, "Load 'curvature' z coord");
        bLoadScalar->setCallback([this] {
            std::cout << "Using z coord as fake curvature...";

            
            // Create curvature function, zero z
            VertexData<double> f(mesh);
            for(VertexPtr v : mesh->vertices()) {
                f[v] = geometry->position(v).z;
                geometry->position(v).z = 0.0;
            }

            // All kinds of things need to be updated
            shapeOpt = new EulerianShapeOptimizer(geometry);
            sceneMeshCuts = new SceneMeshCuts(*(mainGui->viewer), geometry);
            mainGui->replaceMainSceneMesh( sceneMeshCuts );
            mainGui->sceneMesh = sceneMeshCuts;

            // Set
            shapeOpt->setFakeFaceCuvature(f); 

            (*updateVisualization)();
        });
    }

    { // From file
        nanogui::Button *bLoadState = new nanogui::Button(popup, "From file");
        bLoadState->setCallback([this] {
            std::cout << "Selecting state to load...  ";
            
            std::string stateFilename = nanogui::file_dialog({{"msdf", "Multi Signed Distance Function"}}, false);
            if(stateFilename == "") {
                std::cout << "...no file selected" << std::endl;
                return;
            }
            std::cout << "...selected state file " << stateFilename << std::endl;
        
            // Try to load the state. If something goes wrong, catch the exception instead of exiting.
            try {
                shapeOpt->loadFromFile(stateFilename);
                shapeOpt->iIter = 0;
            } catch (const std::exception &exc) {
                std::cout << "State loading failed." << std::endl;
                mainGui->reportError(exc.what());
            } catch (...) {
                std::cout << "State loading failed." << std::endl;
            }
            
            (*updateVisualization)();
        });
    }

    // Add a button to detect symmetry
    nanogui::Button *bSym = new nanogui::Button(window, "Detect symmetry");
    bSym->setCallback([&] {

        SymmetryResult sym = detectSymmetryAuto(geometry);
        // SymmetryResult sym = detectSymmetryAutoRotation(geometry);
        // SymmetryResult sym = detectSymmetryDoubleMirror(geometry);

        if(!sym.symmetryFound) {
            mainGui->reportError("No symmetry found.");
            return;
        }

        shapeOpt->symmetry = sym;
        shapeOpt->hasSymmetry = true;
        shapeOpt->projectSymmetric();
        shapeOpt->clearCaches();

        // Manually swtich the viz to show symmetry        
        vizName = "Symmetry";
        vizPopupPanelButton->setCaption("Symmetry");
        vizPopupPanelButton->setPushed(false);

        // Visualize
        (*updateVisualization)();

    });
    

    new nanogui::Label(window, "Iteration", "sans-bold");

    // Add a button to optimize cuts/region
    nanogui::Button *bStepOpt = new nanogui::Button(window, "Take step");
    bStepOpt->setCallback([&] {

        (*updateParameters)();

        // Do work 
        if(useLineSearchCB->checked()) {
            shapeOpt->doStepLineSearch();
        } else {
            shapeOpt->doStep();
        }


        // Visualize
        (*updateVisualization)();

    });
    
    nanogui::Button *bRunOpt = new nanogui::Button(window, "Take many steps");
    bRunOpt->setCallback([&] {
        
        (*updateParameters)();

        int nSteps = 0;
        tryParseInt(nStepsBox, nSteps);

        // Do work 
        bool movie = makeMovieCB->checked();
        bool rotate = rotateCB->checked();
        bool render = updateRenderCB->checked();
        bool save = saveEveryCB->checked();
        for(int i = 0; i < nSteps; i++) {

            if(useLineSearchCB->checked()) {
                int nSteps = shapeOpt->doStepLineSearch();

                // Tweak the step size param

            } else {
                shapeOpt->doStep();
            }

            if(movie || render) {
                (*updateVisualization)();
                mainGui->viewer->draw();
            }

            if(movie) {
                mainGui->viewer->screenshot();
                // (with ffmpeg: ffmpeg -i screenshot%010d.tga -crf 10 -y movie.mp4)
            }

            // Rotate, if requested
            if(rotate) {
                double delTheta = PI / 200;
                Camera* c = &(mainGui->viewer->camera);
                c->cameraDirection = c->cameraDirection.rotate_around(Vector3{0,1,0}, delTheta);
                c->upDirection = c->upDirection.rotate_around(Vector3{0,1,0}, delTheta);
            }

            // Save, if requested
            if(save) {
            //if(save && shapeOpt->iIter % 100 == 0) {

                { // Save msdf
                    std::string filename = "autosave_" + std::to_string(shapeOpt->iIter) + ".msdf";
                    shapeOpt->saveToFile(filename);
                }


                // { // Save lines
                //     std::string filename = "output_lines_" + std::to_string(shapeOpt->iIter) + ".lines";
                //     shapeOpt->saveCutLines(filename);
                // }

                //{ // Save ply with color
                    //std::string filename = "output_distortion_color_" + std::to_string(shapeOpt->iIter) + ".ply";
                    //saveDistortionPLY(shapeOpt, filename);
                //}

                //{ // Distortion histogram
                    //std::string stateFilename = "distortion_historgram_" + std::to_string(shapeOpt->iIter) + ".txt";
                    //shapeOpt->saveDistortionHistogram(stateFilename);
                //}

                
                // Cylindrical parameterization
                {
                    CornerData<Vector2> paramCoords;
                    Geometry<Euclidean>* cutGeodesicGeom = shapeOpt->getGeodesicFlattenedGeometry(paramCoords); 
    
                    // Process to a disk parameterization
                    CornerData<Vector2> diskCoords;
                    computeCylinderParameterization(cutGeodesicGeom, diskCoords);
    
                    // Save
                    std::string filename = "output_cyl_" + std::to_string(shapeOpt->iIter) + ".obj";
                    WavefrontOBJ::write(filename, *cutGeodesicGeom, diskCoords);    
                }
                
            }

        }

        // Visualize
        (*updateVisualization)();

    });
    
    useLineSearchCB = new nanogui::CheckBox(window, "Line search");

    // Grid for stepping options
	nanogui::Widget* stepOptionWidget = new nanogui::Widget(window);
	nanogui::GridLayout* stepLayout = new nanogui::GridLayout(nanogui::Orientation::Horizontal, 2,
										   nanogui::Alignment::Middle, 15, 5);
	stepLayout->setColAlignment({ nanogui::Alignment::Maximum, nanogui::Alignment::Fill });
	stepLayout->setSpacing(0, 10);
	stepOptionWidget->setLayout(stepLayout);
    
    
    // Stepsize textbox
    new nanogui::Label(stepOptionWidget, "Step size:", "sans-bold");
    stepSizeBox = new nanogui::TextBox(stepOptionWidget);
    initializeTextBox(stepSizeBox, shapeOpt->stepSizeParam);

    
    // nSteps textbox
    new nanogui::Label(stepOptionWidget, "# steps:", "sans-bold");
    nStepsBox = new nanogui::TextBox(stepOptionWidget);
    initializeTextBox(nStepsBox, 10);

    updateRenderCB = new nanogui::CheckBox(window, "Update viz");
    makeMovieCB = new nanogui::CheckBox(window, "Generate movie");
    rotateCB = new nanogui::CheckBox(window, "Auto-rotate");
    saveEveryCB = new nanogui::CheckBox(window, "Save every iteration");



    // Length term
    {
        TermWidget* tWidget = new TermWidget(paramWindow, "Length Regularization", true);
        tWidget->pushParameters = new std::function<void()>([&, tWidget]() {
            tryParseDouble(tWidget->weightBox, shapeOpt->weightLengthRegularization);
            shapeOpt->localScaleLengthRegularization = tWidget->localCB->checked();
        });
        tWidget->weightBox->setValue(std::to_string(1.0));
        termWidgets.push_back(tWidget);
    }
    
    // Bilaplacian term
    {
        TermWidget* tWidget = new TermWidget(paramWindow, "Bilaplacian Regularization", false);
        tWidget->pushParameters = new std::function<void()>([&, tWidget]() {
            tryParseDouble(tWidget->weightBox, shapeOpt->weightBilapRegularization);
        });
        tWidget->weightBox->setValue(std::to_string(0.0));
        termWidgets.push_back(tWidget);
    }
    
    // Dirichlet Distortion term
    {
        TermWidget* tWidget = new TermWidget(paramWindow, "Dirichlet Distortion", true);
        tWidget->pushParameters = new std::function<void()>([&, tWidget]() {
            tryParseDouble(tWidget->weightBox, shapeOpt->weightDirichletDistortion);
            shapeOpt->localScaleDirichletDistortion = tWidget->localCB->checked();
        });
        termWidgets.push_back(tWidget);
    }
    
    // Hencky Distortion term
    {
        TermWidget* tWidget = new TermWidget(paramWindow, "Hencky Distortion");
        tWidget->pushParameters = new std::function<void()>([&, tWidget]() {
            tryParseDouble(tWidget->weightBox, shapeOpt->weightHenckyDistortion);
            // shapeOpt->localScaleHenckyDistortion = tWidget->localCB->checked();
        });
        termWidgets.push_back(tWidget);
    }
    
    // Visibility term
    {
        TermWidget* tWidget = new TermWidget(paramWindow, "Visibility");
        tWidget->pushParameters = new std::function<void()>([&, tWidget]() {
            tryParseDouble(tWidget->weightBox, shapeOpt->weightVisibility);
            // shapeOpt->localScaleVisibility = tWidget->localCB->checked();
        });
        termWidgets.push_back(tWidget);
    }
    
    // Area term
    {
        TermWidget* tWidget = new TermWidget(paramWindow, "Area");
        new nanogui::Label(tWidget->widget, "Region 0 area:", "sans-bold");
        patchZeroAreaBox = new nanogui::TextBox(tWidget->widget);
        initializeTextBox(patchZeroAreaBox, -1.0);
        tWidget->pushParameters = new std::function<void()>([&, tWidget]() {
            tryParseDouble(tWidget->weightBox, shapeOpt->weightArea);
            // shapeOpt->localScaleArea = tWidget->localCB->checked();
            tryParseDouble(patchZeroAreaBox, shapeOpt->patchZeroTargetArea);
        });
        termWidgets.push_back(tWidget);
    }
    
    // Normal deviation term
    {
        TermWidget* tWidget = new TermWidget(paramWindow, "Normal Deviation");
        tWidget->pushParameters = new std::function<void()>([&, tWidget]() {
            tryParseDouble(tWidget->weightBox, shapeOpt->weightNormalDeviation);
            // shapeOpt->localScaleNormalDeviation = tWidget->localCB->checked();
        });
        termWidgets.push_back(tWidget);
    }
    
    

    new nanogui::Label(window, "Visualization", "sans-bold");
    
    
    // Add a radio button to choose what data is being visualized on the mesh
    {
        vizPopupPanelButton = new nanogui::PopupButton(window, vizName);
        nanogui::Popup *popup = vizPopupPanelButton->popup();
        popup->setLayout(new nanogui::GridLayout(nanogui::Orientation::Horizontal, 2,
										         nanogui::Alignment::Fill, 15, 5));
        std::vector<std::string> vizNameOptions = {   "Region Labels",
                                                "Phi",
                                                "Cut Mesh",
                                                "Cut Mesh Param Coarse",
                                                "Cut Mesh Param",
                                                "Distortion",
                                                "Curvature",
                                                "Developable Approximation",
                                                "Patch Area",
                                                "Patch Normal Deviation",
                                                "Boundary Distance",
                                                "Shape Motion",
                                                "Mesh-k",
                                                "Shape Motion-k",
                                                "Feature Size",
                                                "Scaled Distortion",
                                                "Visibility",
                                                "Symmetry",
                                                "DEBUG",
                                                "None"};
        for(std::string name : vizNameOptions) {
            nanogui::Button *b = new nanogui::Button(popup, name);
            b->setFixedWidth(200);
            b->setCallback([&, name] {
                vizName = name;
                vizPopupPanelButton->setCaption(name);
                vizPopupPanelButton->setPushed(false);
                (*updateVisualization)();
            });
        }
    }

    // new nanogui::Label(window, "Vector Visualization", "sans-bold");


    // Add a radio button to choose what data is being visualized as vectors
    /*
    vectorVisualizeCombo = new nanogui::ComboBox(window, {"du/dn", "du/dn Smoothed", "Shape Optimizer", "Regularizer", "Update", "None"});
    vectorVisualizeCombo->setCallback([this](int newState) {
        (*updateVisualization)();
    });
    vectorVisualizeCombo->setSelectedIndex(vectorVisualizeCombo->items().size()-1); // start on "None"
    */

    // Grid for viz options
	nanogui::Widget* vizOptionWidget= new nanogui::Widget(window);
	nanogui::GridLayout* vizLayout = new nanogui::GridLayout(nanogui::Orientation::Horizontal, 2,
										   nanogui::Alignment::Middle, 15, 5);
	vizLayout->setColAlignment({ nanogui::Alignment::Maximum, nanogui::Alignment::Fill });
	vizLayout->setSpacing(0, 10);
	vizOptionWidget->setLayout(vizLayout);
    
    // Label to show how many patches there are
	new nanogui::Label(vizOptionWidget, "# Patches:", "sans-bold");
	nPatchLabel = new nanogui::Label(vizOptionWidget, "-1", "sans-bold");

    // Textbox for picking a region to show data for
	new nanogui::Label(vizOptionWidget, "Inspect patch/region:", "sans-bold");
	kInspectBox = new nanogui::TextBox(vizOptionWidget);
    initializeTextBox(kInspectBox , 0);
    kInspectBox->setCallback([this](const std::string &newVal) {
        (*updateVisualization)();
        return true;
    });

    // Add a checkbox to show boundary lines
    showLinesCB = new nanogui::CheckBox(window, "Show boundary");
    showLinesCB->setCallback([this](bool newState) {
        (*updateVisualization)();
    });
    
    
    // Add a checkbox to show extra lines
    showExtraLinesCB = new nanogui::CheckBox(window, "Show extra cuts");
    showExtraLinesCB->setCallback([this](bool newState) {
        (*updateVisualization)();
    });
    

    // === Input/Output
    new nanogui::Label(window, "Input/Output", "sans-bold");



    // Popup to save various things
    nanogui::PopupButton *popupSaveButton = new nanogui::PopupButton(window, "Save", ENTYPO_ICON_EXPORT);
    popup = popupSaveButton->popup();
    popup->setLayout(new nanogui::GroupLayout());

    { // Distortion histogram
        nanogui::Button *bSaveHist = new nanogui::Button(popup, "Save Distortion Histogram");
        bSaveHist->setCallback([this] {
            std::string stateFilename = "distortion_historgram_" + std::to_string(shapeOpt->iIter) + ".txt";
            shapeOpt->saveDistortionHistogram(stateFilename);
        });
    }

    { // Save msdf
        nanogui::Button *bSaveState = new nanogui::Button(popup, "Save MSDF");
        bSaveState->setCallback([this] {
            std::cout << "Selecting state to save...  ";
            
            std::string stateFilename = nanogui::file_dialog({{"msdf", "Multi Signed Distance Function"}}, true);
            if(stateFilename == "") {
                std::cout << "...no file selected" << std::endl;
                return;
            }
            std::cout << "...selected state file " << stateFilename << std::endl;
        
            // Try to save the state. If something goes wrong, catch the exception instead of exiting.
            try {
                shapeOpt->saveToFile(stateFilename);
            } catch (const std::exception &exc) {
                std::cout << "State saving failed." << std::endl;
                mainGui->reportError(exc.what());
            } catch (...) {
                std::cout << "State saving failed." << std::endl;
            }
        });
    }

    
    { // Save obj file with coordinates
        nanogui::Button *bSaveState = new nanogui::Button(popup, "Save .obj with injective texcoords");
        bSaveState->setCallback([this] {
            std::cout << "Selecting state to save...  ";
            
            std::string stateFilename = nanogui::file_dialog({{"obj", "Wavefront OBJ"}}, true);
            if(stateFilename == "") {
                std::cout << "...no file selected" << std::endl;
                return;
            }
            std::cout << "...selected state file " << stateFilename << std::endl;
        
            // Try to save the state. If something goes wrong, catch the exception instead of exiting.
            try {

                // Get the mesh and parameterization
                CornerData<Vector2> paramCoords;
                Geometry<Euclidean>* cutGeodesicGeom = shapeOpt->getGeodesicFlattenedGeometry(paramCoords); 

                // Get region indices of faces
                FaceData<int> faceComponents = computeFaceComponents(cutGeodesicGeom, paramCoords);

                // Save
                // WavefrontOBJ::write(stateFilename, *cutGeodesicGeom, paramCoords, faceComponents);
                WavefrontOBJ::write(stateFilename, *cutGeodesicGeom, paramCoords);

            } catch (const std::exception &exc) {
                std::cout << "State saving failed." << std::endl;
                mainGui->reportError(exc.what());
            } catch (...) {
                std::cout << "State saving failed." << std::endl;
            }
        });
    }
    
    { // Save boundary as svg
        nanogui::Button *bSaveState = new nanogui::Button(popup, "Save boundary as svg");
        bSaveState->setCallback([this] {
            std::cout << "Selecting state to save...  ";
            
            std::string stateFilename = nanogui::file_dialog({{"svg", "Silly Vector Graphics"}}, true);
            if(stateFilename == "") {
                std::cout << "...no file selected" << std::endl;
                return;
            }
            std::cout << "...selected state file " << stateFilename << std::endl;
        
            // Try to save the state. If something goes wrong, catch the exception instead of exiting.
            try {

                // Get the mesh and parameterization
                CornerData<Vector2> paramCoords;
                Geometry<Euclidean>* cutGeodesicGeom = shapeOpt->getGeodesicFlattenedGeometry(paramCoords); 

                // Save
                writeBoundarySVG(stateFilename, cutGeodesicGeom, paramCoords);    

            } catch (const std::exception &exc) {
                std::cout << "State saving failed." << std::endl;
                mainGui->reportError(exc.what());
            } catch (...) {
                std::cout << "State saving failed." << std::endl;
            }
        });
    }

    { // Save obj file with coordinates
        nanogui::Button *bSaveState = new nanogui::Button(popup, "Save .obj with texcoord as position");
        bSaveState->setCallback([this] {
            std::cout << "Selecting state to save...  ";
            
            std::string stateFilename = nanogui::file_dialog({{"obj", "Wavefront OBJ"}}, true);
            if(stateFilename == "") {
                std::cout << "...no file selected" << std::endl;
                return;
            }
            std::cout << "...selected state file " << stateFilename << std::endl;
        
            // Try to save the state. If something goes wrong, catch the exception instead of exiting.
            try {

                // Get the mesh and parameterization
                CornerData<Vector2> paramCoords;
                Geometry<Euclidean>* cutGeodesicGeom = shapeOpt->getGeodesicFlattenedUVPosGeometry(paramCoords); 

                // Save
                WavefrontOBJ::write(stateFilename, *cutGeodesicGeom, paramCoords);    

            } catch (const std::exception &exc) {
                std::cout << "State saving failed." << std::endl;
                mainGui->reportError(exc.what());
            } catch (...) {
                std::cout << "State saving failed." << std::endl;
            }
        });
    }
    
    { // Save extrinsic developable obj
        nanogui::Button *bSaveState = new nanogui::Button(popup, "Save .obj extrinsic developable");
        bSaveState->setCallback([this] {
            std::cout << "Selecting state to save...  ";
            
            std::string stateFilename = nanogui::file_dialog({{"obj", "Wavefront OBJ"}}, true);
            if(stateFilename == "") {
                std::cout << "...no file selected" << std::endl;
                return;
            }
            std::cout << "...selected state file " << stateFilename << std::endl;
        
            // Try to save the state. If something goes wrong, catch the exception instead of exiting.
            try {

                // Get the mesh 
                Geometry<Euclidean>* cutGeodesicGeom = shapeOpt->getExtrinsicDevelopableGeometry();
                
                // Disk parameterize 
                CornerData<Vector2> diskCoords;
                computeCylinderParameterization(cutGeodesicGeom, diskCoords);

                // Save
                WavefrontOBJ::write(stateFilename, *cutGeodesicGeom, diskCoords);    

            } catch (const std::exception &exc) {
                std::cout << "State saving failed." << std::endl;
                mainGui->reportError(exc.what());
            } catch (...) {
                std::cout << "State saving failed." << std::endl;
            }
        });
    }
    
    { // Save disk-parameterized (r-theta) obj
        nanogui::Button *bSaveState = new nanogui::Button(popup, "Save .obj cylinder parameterized");
        bSaveState->setCallback([this] {
            std::cout << "Selecting state to save...  ";
            
            std::string stateFilename = nanogui::file_dialog({{"obj", "Wavefront OBJ"}}, true);
            if(stateFilename == "") {
                std::cout << "...no file selected" << std::endl;
                return;
            }
            std::cout << "...selected state file " << stateFilename << std::endl;
        
            // Try to save the state. If something goes wrong, catch the exception instead of exiting.
            try {

                // Get the mesh and parameterization
                CornerData<Vector2> paramCoords;
                Geometry<Euclidean>* cutGeodesicGeom = shapeOpt->getGeodesicFlattenedGeometry(paramCoords); 

                // Process to a disk parameterization
                CornerData<Vector2> diskCoords;
                computeCylinderParameterization(cutGeodesicGeom, diskCoords);

                // Save
                WavefrontOBJ::write(stateFilename, *cutGeodesicGeom, diskCoords);    

            } catch (const std::exception &exc) {
                std::cout << "State saving failed." << std::endl;
                mainGui->reportError(exc.what());
            } catch (...) {
                std::cout << "State saving failed." << std::endl;
            }
        });
    }
    
    { // Save  lines
        nanogui::Button *bSaveState = new nanogui::Button(popup, "Save cut lines");
        bSaveState->setCallback([this] {
            std::cout << "Selecting state to save...  ";
            
            std::string stateFilename = nanogui::file_dialog({{"lines", "Lines (ASCII)"}}, true);
            if(stateFilename == "") {
                std::cout << "...no file selected" << std::endl;
                return;
            }
            std::cout << "...selected state file " << stateFilename << std::endl;
        
            // Try to save the state. If something goes wrong, catch the exception instead of exiting.
            try {
                shapeOpt->saveCutLines(stateFilename);
            } catch (const std::exception &exc) {
                std::cout << "State saving failed." << std::endl;
                mainGui->reportError(exc.what());
            } catch (...) {
                std::cout << "State saving failed." << std::endl;
            }
        });
    }

    { // Save ply with colors
        nanogui::Button *bSaveState = new nanogui::Button(popup, "Save .ply with distortion color");
        bSaveState->setCallback([this] {
            std::cout << "Selecting file to save...  ";
            
            std::string saveFilename = nanogui::file_dialog({{"ply", "PLY"}}, true);
            if(saveFilename == "") {
                std::cout << "...no file selected" << std::endl;
                return;
            }
            std::cout << "...selected file " << saveFilename << std::endl;
        
            // Try to save the state. If something goes wrong, catch the exception instead of exiting.
            try {
                saveDistortionPLY(shapeOpt, saveFilename);

            } catch (const std::exception &exc) {
                std::cout << "State saving failed." << std::endl;
                mainGui->reportError(exc.what());
            } catch (...) {
                std::cout << "State saving failed." << std::endl;
            }
        });
    }

   mainGui->viewer->mainScreen->performLayout();

   // Initialize visualiation
   (*updateVisualization)();


   return true; // success
}


void CutsGui::hide(void) {
    delete sceneMeshCuts; sceneMeshCuts = nullptr;
    delete shapeOpt; shapeOpt = nullptr;
    delete updateVisualization; updateVisualization = nullptr;
    delete updateParameters; updateParameters = nullptr;

    paramWindow->dispose();
    paramWindow = nullptr;

    ProjectGui::hide();
}

void CutsGui::revert(void)
{}

void CutsGui::updatePickedElement(const int pickedID, int x, int y) {

    cout << endl << "Pick called with " << pickedID << endl;

    /*
    if(!mesh) return;

    size_t nEdge = mesh->nEdges();
    size_t nFace = mesh->nFaces();
    size_t nVert = mesh->nVertices();

    if(pickedID < (int)(nEdge+nFace)) return;
    if(pickedID >= (int)(nEdge+nFace+nVert)) return;

    size_t vID = pickedID - nEdge - nFace;
    VertexPtr v = mesh->vertex(vID);
    
    cout << endl << "Picked vertex " << v << endl;
    cout << "phi[v] = " << endl;
    cout << shapeOpt->phi[v] << endl;
            
    // cout << "mean = " << mean << " smallest = " << smallestVal << endl;
    cout << "zFormDiag = " << shapeOpt->zeroFormLaplacian(shapeOpt->vInd[v], shapeOpt->vInd[v]) << endl;
    cout << "regularizer = " << (shapeOpt->laplaceRegularizer[v]) << endl;

    cout << "Neighbors: " << endl;
    for(VertexPtr vN : v.adjacentVertices()) {
        cout << endl << vN << endl << shapeOpt->phi[vN] << endl;
        cout << "zFormE = " << shapeOpt->zeroFormLaplacian(shapeOpt->vInd[v], shapeOpt->vInd[vN]) << endl;
    }
    */
    
}

