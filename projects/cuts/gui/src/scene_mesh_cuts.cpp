#include "scene_mesh_cuts.h"

#include "halfedge_mesh.h"
#include "viewer.h"
#include "direction_fields.h"

#include "colormaps.h"
#include "shaders.h"
#include "shaders/shiny_shaders.h"
#include "shaders/checker_shaders.h"
#include "cut_shaders.h"

// Only so that we can customize the region visualization for the number of regions
#include "eulerian_shape_optimizer.h"

SceneMeshCuts::SceneMeshCuts(Viewer& parent_, Geometry<Euclidean>* geometryToDraw) :
    SceneMesh(parent_, geometryToDraw),
    vectorArtist(parent_), lineArtist(parent_), lineArtistAlternate(parent_)
{

    shapeLengthScale = geometry->lengthScale();

    vectorArtist.enabled = false;
    lineArtist.enabled = false;
    lineArtistAlternate.enabled = false;

    lineArtist.lineRadius = 0.001 * shapeLengthScale;
    lineArtist.lineColor = RGB_RED;
    lineArtistAlternate.lineRadius = 0.001 * shapeLengthScale;
    lineArtistAlternate.lineColor = RGB_ORANGE;
}
        
SceneMeshCuts::~SceneMeshCuts() {
    showNothing();
}

void SceneMeshCuts::draw() {

    vectorArtist.draw();

    SceneMesh::draw();
}

void SceneMeshCuts::setCurrentMesh(Geometry<Euclidean>* geometry_) {
    geometry = geometry_;
    mesh = geometry->getMesh();
    fillMeshBuffers();
    fillMeshBuffersShadow();
    fillMeshBuffersPick();
}


Vector3 SceneMeshCuts::getPaletteColor(size_t iRegion) {
    return CM_SPECTRAL.getValue(((double)iRegion) / (K_REGIONS-1));
}

void SceneMeshCuts::showScalar(const VertexData<double>& vals, bool symmetric) {
    showNothing();
    
    if(symmetric) {
        double min, max;
        computeSymmetricBounds(vals, min, max);
        setColorFromScalar(vals, min, max, CM_COOLWARM);
    } else {
        setColorFromScalar(vals);
    }

}

void SceneMeshCuts::showScalar(const FaceData<double>& vals, bool symmetric) {
    showNothing();
    
    if(symmetric) {
        double min, max;
        computeSymmetricBounds(vals, min, max);
        setColorFromScalar(vals, min, max, CM_COOLWARM);
    } else {
        setColorFromScalar(vals);
    }

}

void SceneMeshCuts::showSymmetry(const SymmetryResult& s) {
    showNothing();

    if(!s.symmetryFound) {
        return;
    }

    Vector3 amber = Vector3{255, 193, 7}/255;
    Vector3 blue = Vector3{33, 150, 243}/255;

    VertexData<Vector3> isCanonicalVert(mesh, blue);
    for(VertexPtr v : s.canonicalVertices) {
        isCanonicalVert[v] = amber;
    }
    setColorFromRGB(isCanonicalVert);
}

// Visualize all-positive linear combinations of values
void SceneMeshCuts::showMultiScalar(const FaceData<Eigen::VectorXd>& vals) {

    // First compute a scale, which is the largest sum
    double maxSum = 0;
    for(FacePtr f : mesh->faces()) {
        double sum = vals[f].cwiseAbs().sum();
        maxSum = std::max(maxSum, sum);
    }

    // Each color is a linear combination of the palette
    FaceData<Vector3> resultColor(mesh, Vector3{0.0, 0.0, 0.0});
    for(FacePtr f : mesh->faces()) {
        for(size_t i = 0; i < (size_t)vals[f].rows(); i++) {
            double ratio = clamp(vals[f](i), 0.0, std::numeric_limits<double>::infinity()) / maxSum;
            resultColor[f] += ratio * getPaletteColor(i);
        }

        // Fill with white (rather than the implicit black)
        double mySum = vals[f].cwiseAbs().sum();
        resultColor[f] += (maxSum - mySum) / maxSum * Vector3{1.0, 1.0, 1.0};
    }
    
    setColorFromRGB(resultColor);
}

// Visualize all-positive linear combinations of values
void SceneMeshCuts::showMultiScalar(const VertexData<Eigen::VectorXd>& vals) {

    // First compute a scale, which is the largest sum
    double maxSum = 0;
    double meanSum = 0;
    for(VertexPtr v : mesh->vertices()) {
        double sum = vals[v].cwiseAbs().sum();
        meanSum += sum;
        maxSum = std::max(maxSum, sum);
    }
    meanSum /= mesh->nVertices();
    cout << "Max sum = " << maxSum << endl;

    // Each color is a linear combination of the palette
    VertexData<Vector3> resultColor(mesh, Vector3{0.0, 0.0, 0.0});
    for(VertexPtr v : mesh->vertices()) {
        double mySum = vals[v].cwiseAbs().sum();
        double mag = std::max(mySum,meanSum);
        for(size_t i = 0; i < (size_t)vals[v].rows(); i++) {
            // double ratio = clamp(vals[v](i), 0.0, std::numeric_limits<double>::infinity()) / maxSum; FIXME
            double ratio = clamp(vals[v](i), 0.0, std::numeric_limits<double>::infinity()) / mag;
            resultColor[v] += ratio * getPaletteColor(i);
        }

        // Fill with white (rather than the implicit black)
        resultColor[v] += (mag - mySum) / mag * Vector3{1.0, 1.0, 1.0};
    }
    
    setColorFromRGB(resultColor);
}

void SceneMeshCuts::showRegionLabels(const FaceData<int>& regionLabels) {
    showNothing();

    cout << "Showing regions" << endl;
    

    // Map integer labels to colors
    FaceData<Vector3> faceRGB(mesh);
    for(FacePtr f : mesh->faces()) {
        faceRGB[f] = getPaletteColor(regionLabels[f]);
    }

    setColorFromRGB(faceRGB);
}

void SceneMeshCuts::showRegionLabels(const VertexData<Eigen::VectorXd>& vals) {
    showNothing();

    cout << "Showing regions" << endl;

    delete faceGLProgram;
    faceGLProgram = new GLProgram(&MULTI_REGION_VERT_SHADER, &MULTI_REGION_FRAG_SHADER, DrawMode::Triangles);


    std::vector<Vector3> regionVals(3*mesh->nFaces());
    std::array<std::vector<Vector3>,3> regionColors = {{std::vector<Vector3>(3*mesh->nFaces()),
                                                        std::vector<Vector3>(3*mesh->nFaces()),
                                                        std::vector<Vector3>(3*mesh->nFaces())}};

    size_t i = 0;
    for(FacePtr f : mesh->faces()) {

        std::array<VertexPtr,3> vNeigh = {{ f.halfedge().vertex(),
                                            f.halfedge().next().vertex(),
                                            f.halfedge().next().next().vertex()}};

        std::array<int,3> minLabel;
        for(size_t j = 0; j < 3; j++) {
            int myRegion;
            vals[vNeigh[j]].minCoeff(&myRegion);
            minLabel[j] = myRegion;
        }

        for(size_t j = 0; j < 3; j++) {

            for(size_t k = 0; k < 3; k++) {
                regionColors[k][i] = getPaletteColor(minLabel[k]);
                regionVals[i][k] = vals[vNeigh[j]][minLabel[k]];
            }

            i++;
        }

        

    }

    faceGLProgram->setAttribute("a_regionVals", regionVals);
    faceGLProgram->setAttribute("a_regionColor0", regionColors[0]);
    faceGLProgram->setAttribute("a_regionColor1", regionColors[1]);
    faceGLProgram->setAttribute("a_regionColor2", regionColors[2]);

    fillMeshBuffers();
}
void SceneMeshCuts::showRegionChecker(FaceData<int> regionLabels, CornerData<Vector2> params) {
    showNothing();

    cout << "Showing checker regions" << endl;

    delete faceGLProgram;
    faceGLProgram = new GLProgram(&CUT_CHECKER_VERT_SHADER, &CUT_CHECKER_FRAG_SHADER, DrawMode::Triangles);

    std::vector<Vector2> UVs(3 * mesh->nFaces());
    std::vector<Vector3> colors(3 * mesh->nFaces());

    size_t i = 0;
    for(FacePtr f : mesh->faces()) {
        for(HalfedgePtr he : f.adjacentHalfedges()) {
            colors[i] = getPaletteColor(regionLabels[f]);
            UVs[i] = params[he.next().corner()] / 25;
            i++;
        }
    }

    faceGLProgram->setAttribute("a_color", colors);
    faceGLProgram->setAttribute("a_texcoord", UVs);

    // faceGLProgram->setUniform("u_intensity", 1.0);

    fillMeshBuffers();
}

        
void SceneMeshCuts::showRegionChecker(CornerData<Vector2> params) {
    showNothing();

    cout << "Showing checker regions" << endl;

    delete faceGLProgram;
    faceGLProgram = new GLProgram(&CUT_CHECKER_VERT_SHADER, &CUT_CHECKER_FRAG_SHADER, DrawMode::Triangles);

    std::vector<Vector2> UVs(3 * mesh->nFaces());
    std::vector<Vector3> colors(3 * mesh->nFaces());

    size_t i = 0;
    for(FacePtr f : mesh->faces()) {
        for(HalfedgePtr he : f.adjacentHalfedges()) {
            colors[i] = RGB_LIGHTGRAY;
            UVs[i] = params[he.next().corner()] / (shapeLengthScale * 5);
            i++;
        }
    }

    faceGLProgram->setAttribute("a_color", colors);
    faceGLProgram->setAttribute("a_texcoord", UVs);

    // faceGLProgram->setUniform("u_intensity", 1.0);

    fillMeshBuffers();
}

void SceneMeshCuts::setLines(const std::vector<std::array<Vector3,2>>& lines) {
    lineArtist.setLines(lines);
}

void SceneMeshCuts::toggleLines(bool b) {
    lineArtist.enabled = b;
}

void SceneMeshCuts::setLinesAlternate(const std::vector<std::array<Vector3,2>>& lines) {
    lineArtistAlternate.setLines(lines);
}

void SceneMeshCuts::toggleLinesAlternate(bool b) {
    lineArtistAlternate.enabled = b;
}

void SceneMeshCuts::showNothing() {

    // Clear out the vector artist
    vectorArtist.enabled = false;

    setColorDefault();
}

void SceneMeshCuts::showVector(const FaceData<Vector3>& vals) {
    showNothingVector();


    vectorArtist.setVectors(geometry, vals, true);
    // std::vector<Vector3> fakeV = {Vector3{1,1,1}};
    // vectorArtist.setVectors(fakeV,fakeV);
    vectorArtist.enabled = true;
}

void SceneMeshCuts::showOneForm(const EdgeData<double>& vals, bool sparse) {
    showNothingVector();

    vectorArtist.setOneForm(geometry, vals, true, sparse);
    vectorArtist.enabled = true;
}

void SceneMeshCuts::showEdgePerpVector(const EdgeData<double>& vals) {
    showNothingVector();

    vectorArtist.setEdgePerpVectors(geometry, vals, true);
    vectorArtist.enabled = true;
    
}

void SceneMeshCuts::showEdgeVector(const EdgeData<double>& vals) {
    showNothingVector();

    vectorArtist.setEdgeVectors(geometry, vals, true);
    vectorArtist.enabled = true;
    
}

void SceneMeshCuts::showNothingVector() {
    // Clear out the vector artist
    vectorArtist.enabled = false;

}

void SceneMeshCuts::computeSymmetricBounds(VertexData<double> data, double& min, double &max, double cap) {

    double extreme = 0;
    for(VertexPtr v : mesh->vertices()) {
        extreme = std::max(extreme, std::abs(data[v]));
    }

    extreme = std::min(extreme, cap);

    if(extreme == 0) {
        extreme = 0.001;
    }

    min = -extreme;
    max = extreme;

    // cout << "min = " << min << " max = " << max << endl;
}

void SceneMeshCuts::computeSymmetricBounds(FaceData<double> data, double& min, double &max, double cap) {

    double extreme = 0;
    for(FacePtr f : mesh->faces()) {
        extreme = std::max(extreme, std::abs(data[f]));
    }

    extreme = std::min(extreme, cap);

    if(extreme == 0) {
        extreme = 0.001;
    }

    min = -extreme;
    max = extreme;

    // cout << "min = " << min << " max = " << max << endl;
}