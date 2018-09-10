#pragma once

#include <vector>
#include <cmath>
#include <array>
#include <unordered_map>

#include "utilities.h"
#include "geometry.h"
#include "surface_patch.h"
#include "fast_cholesky.h"
#include "fast_triangle_soup.h"
#include "disjoint_sets.h"
#include "detect_symmetry.h"

#define K_REGIONS 10

typedef Eigen::Matrix<double, K_REGIONS, 1> LabelVec; 


// Utility class for representing a region boundary embedded in the mesh
enum class BType {TRI, QUAD, TRIPLE};
struct BoundarySegment {

    BType type;
    HalfedgePtr heStart, heEnd; // following he.next() traces the region bounded by this segment
    double tStart, tEnd;        // parameterized along the halfedges above
    // size_t index;               // used for building matrices

    Vector3 triplePoint;         // for type=TRIPLE, the point where the three regions meet
    // Vector3 tripleBary;          // barycentric coordinates of the triple point, starting at face.he().vertex()

};

class EulerianShapeOptimizer {
    
public:

    EulerianShapeOptimizer(Geometry<Euclidean> *geometry);
    ~EulerianShapeOptimizer();


    // === State ===

    // Primary state
    int iIter = 0;
    VertexData<LabelVec> phi; // signed distance function --> negative for current region, positive for others

    // Parameters
    double stepSizeParam = 1.0; // note that this is just saved state, actual step is passed in to takeGradientStep()

    double weightLengthRegularization = 1.0;
    double weightBilapRegularization = 1.0;
    double weightDirichletDistortion = 0.0;
    double weightHenckyDistortion = 0.0;
    double weightVisibility = 0.0;
    double weightArea = 0.0;
    double weightNormalDeviation = 0.0;
    
    bool localScaleLengthRegularization = false;
    bool localScaleDirichletDistortion = false;
    // bool localScaleHenckyDistortion = false;
    bool localScaleVisibility = false;
    // bool localScaleArea = false;
    // bool localScaleNormalDeviation = false;
    double patchZeroTargetArea = -1.0;

    // Symmetry information
    bool hasSymmetry = false;
    SymmetryResult symmetry;
    void projectSymmetric();
    
    // = Space filling curve stuff 
    // don't forget to check if we are using constant curvature or not!
    bool spaceFillingCurveMode = false;
    double stepSizeFactor = 1.;

    // = Pinned patch stuff
    bool pinPatchMode = false;
    size_t pinVert = 32958; // what vertex do we require to lie in the patch?
    // double minBoundaryDist = .2; // vertex wants to be at least this far from boundary
    double targetBoundaryDist = .1; // vertex wants to be at least this far from boundary
    VertexData<double> generateLocationTerm(); // generate a scalar function (we abuse "visibility" in code names) that measures nearness to the target point


    // === WELCOME TO THE DATA ZOO ===

    // Depending on the parameters used for the optimization, some subset of the values are actually necessary. It is prohibitively
    // expensive to compute them all on every iteration, so instead we use a lazy approach.
    // Each quantity is governed by the corresponding ensureHave___() method. Call this to populate all of the relevant fields listed
    // below it. The ensureHave___() method will only recompute the values if needed (AKA if phi has changed since the last time they
    // were computed). All necessary dependencies are managed internally, one ensureHave___() method will call others as necessary.

    // == Core boundary geometry
    void ensureHaveBoundaryGeometry();
    bool haveBoundaryGeometry = false;
    VertexData<int> regionLabels;
    std::vector<std::vector<BoundarySegment>> boundarySegments;
    HalfedgeData<BoundarySegment*> crossingInnerSegment;   // boundary segment that crosses this halfedge (for which he.vertex() is inside the region)
    HalfedgeData<BoundarySegment*> crossingOuterSegment;   // boundary segment that crosses this halfedge (for which he.vertex() is outside the region)
    FaceData<std::array<BoundarySegment*, 3>> faceTripleSegments;   
    EdgeData<char> boundaryEdgeCrossing; // does a boundary cross this edge?
    FaceData<char> boundaryFaceCrossing; // does a boundary cross through this face?
    void inferVertexRegions();
    void inferBoundaryGeometry();
    void recomputeSDFFromBoundary();
    void snapSDF();

    // == Connected components
    void ensureHaveConnectedComponents();
    bool haveConnectedComponents = false;
    VertexData<int> connectedComponent;
    std::vector<int> componentRegion;
    size_t nConnectedComponents;
    void findConnectedComponents();


    // == Triangle soup cut meshes (main mesh and patch meshes)
    void ensureHaveTriangleSoupCutMeshes();
    bool haveTriangleSoupCutMeshes = false;

    // Cut mesh 
    FastTriangleSoup* cutMesh = nullptr;
    std::vector<TriBool> cutMeshIsComponentBoundaryEdge;
    VertexData<size_t> origMeshToCutMesh;
    std::vector<VertexPtr> cutMeshVertexParent;
    std::vector<size_t> cutMeshFacePatchIndex;
    std::vector<LabelVec> cutMeshInterpolatedPhi;
    std::vector<bool> cutMeshRegionBoundary;
    void extractCutMeshes();

    // Patch meshes
    std::vector<SurfacePatch*> patches;
    void buildPatchMeshes();
    void copyFeatureSizeToPatchMeshes();


    // == Remeshed patches
    // void ensureHaveRemeshedPatches();
    // bool haveRemeshedPatches = false;
    // void remeshPatches();


    // == Poisson problem for distortion shape optimization
    void ensureHaveYamabeSolution();
    bool haveYamabeSolution = false;
    std::vector<double> distortion;
    std::vector<LabelVec> yamabeDuDn;
    std::array<std::vector<double>, K_REGIONS> regionDistortion;
    void solveYamabeProblems();


    // == Visibility term
    void ensureHaveVisibilityFraction();
    bool haveVisibilityFraction = false;
    bool visibilityCasted = false;
    VertexData<double> visibilityFraction;
    std::vector<double> cutMeshVisibility;


    // == Area term
    void ensureHavePatchArea();
    bool havePatchArea = false;
    std::vector<double> connectedComponentArea;
    void computePatchArea();


    // == Normals term
    void ensureHavePatchNormals();
    bool havePatchNormals = false;
    std::vector<Vector3> connectedComponentNormal;
    void computePatchNormals();


    // == Cut halfedge mesh
    void ensureHaveCutHalfedgeMesh();
    bool haveCutHalfedgeMesh = false;
    Geometry<Euclidean>* cutMeshGeometry = nullptr;
    void constructCutHalfedgeMesh();


    // == Halfedge mesh of each patch 
    void ensureHavePatchGeometries();
    bool havePatchGeometries = false;
    void constructPatchGeometries();


    // == Coarse parameterized mesh 
    // (This version has non-disk-like regions cut using edge paths)
    // The relevant mesh is cutMeshGeometry
    void ensureHaveCoarseParameterizedMesh();
    bool haveCoarseParameterizedMesh = false;
    CornerData<Vector2> coarseGeometryCoords;
    void computeCoarseFlattenedGeometry();
    
    
    // == Geodesic parameterized mesh 
    // (This version has s cut using geodesic lines)
    // The relevant mesh is cutMeshGeometry
    void ensureHaveGeodesicParameterizedMesh();
    bool haveGeodesicParameterizedMesh = false;
    Geometry<Euclidean>* geodesicCutMeshGeometry = nullptr;
    CornerData<Vector2> geodesicGeometryCoords;
    void computeGeodesicFlattenedGeometry();
    
    // == Geodesc parameterized mesh with UV in postion
    // (This version has s cut using geodesic lines)
    // The relevant mesh is cutMeshGeometry
    void ensureHaveGeodesicParameterizedUVPosMesh();
    bool haveGeodesicParameterizedUVPosMesh = false;
    Geometry<Euclidean>* geodesicCutMeshUVPosGeometry = nullptr;
    CornerData<Vector2> geodesicGeometryCoordsUVPos;
    void computeGeodesicFlattenedGeometryUVPos();
    
    // == Extrinsic developable mesh
    // (This version has s cut using geodesic lines)
    // The relevant mesh is cutMeshGeometry
    void ensureHaveExtrinsicDevelopableMesh();
    bool haveExtrinsicDevelopableMesh = false;
    Geometry<Euclidean>* extrinsicDevelopableMeshGeometry = nullptr;
    void computeExtrinsicDevelopableMesh();

    
    // === Build the gradient and step along it ===
    void initializeBoundaryGradient();
    void buildGradientAtBoundary(); // builds the explicit terms only (not the boundary length terms)
    void addBoundaryGradientTermLengthRegularization();
    void addBoundaryGradientTermDirichletDistortion();
    void addBoundaryGradientTermHenckyDistortion();
    void addBoundaryGradientTermArea();
    void addBoundaryGradientTermVisibility();
    void addBoundaryGradientTermNormalDeviation();

    
    VertexData<LabelVec> extendGradientToSurface();
    double computeGradientSlope();
    void takeGradientStep(VertexData<LabelVec> gradient, double deltaTParam);
    std::vector<LabelVec> extendedGradient; // save for viz

    // === Compute energy
    double computeEnergy(bool print=false);
    double evaluateEnergyTermLengthRegularization();
    double evaluateEnergyTermDirichletDistortion();
    double evaluateEnergyTermHenckyDistortion();
    double evaluateEnergyTermArea();
    double evaluateEnergyTermVisibility();
    double evaluateEnergyTermNormalDeviation();

    // === Other Methods ===

    // Initialization
    void initializeData();
    void initializeTimestep();
    void computeFeatureSize();
    void setFakeFaceCuvature(VertexData<double> scalarFunction);

    // Iterative optimization
    void clearCaches();
    // bool closeHoles();

    // High level control
    void doStep();
    int doStepLineSearch();
    void solveSystems();
    void setState(VertexData<LabelVec> newPhi);

    // === Input/Output
    void saveToFile(std::string filename);
    void loadFromFile(std::string filename);
    void saveDistortionHistogram(std::string filename);
    void saveCutLines(std::string filename);


    // DEBUG
    // Temporary for debugging
    VertexData<double> getDebug();
    std::vector<double> debugVal;
    VertexData<double> debugValVert;


    // === Visualization Helpers
    VertexData<Eigen::VectorXd> getRegionValues();
    FaceData<double> getPatchArea();
    FaceData<double> getPatchNormalDeviation();
    VertexData<double> getDistortion();
    VertexData<double> getScaledDistortion();
    VertexData<double> getBoundaryDistance();
    VertexData<double> getVisibility();
    Geometry<Euclidean>* getCustomGeometry();
    Geometry<Euclidean>* getCoarseFlattenedGeometry(CornerData<Vector2>& paramCoords);
    Geometry<Euclidean>* getGeodesicFlattenedGeometry(CornerData<Vector2>& paramCoords);
    Geometry<Euclidean>* getGeodesicFlattenedUVPosGeometry(CornerData<Vector2>& paramCoords);
    Geometry<Euclidean>* getExtrinsicDevelopableGeometry();
    Geometry<Euclidean>* getPatchGeometry(int iPatch); 
    VertexData<double> getShapeOptimizerMotion();
    VertexData<double> getRegionShapeOptimizerMotion(int iRegion);
    VertexData<double> getPhiK(int iRegion);
    void getConnectedComponentMeshes(std::vector<Geometry<Euclidean>*>& geometries);
    VertexData<Eigen::VectorXd> getUpdate();
    std::vector<std::array<Vector3, 2>> getBoundaryLines();
    std::vector<std::array<Vector3, 2>> getExtraCutLines();
    Geometry<Euclidean>* generateParameterizedNiceGeometry(CornerData<Vector2>& paramCoords, bool splitPatchSeams=false);

    
    // Hide copy and move constructors, we don't wanna mess with that
    EulerianShapeOptimizer(const EulerianShapeOptimizer &other) = delete;
    EulerianShapeOptimizer &operator=(const EulerianShapeOptimizer &other) = delete;
    EulerianShapeOptimizer(EulerianShapeOptimizer &&other) = delete;
    EulerianShapeOptimizer &operator=(EulerianShapeOptimizer &&other) = delete;

// private:
    // === Private data ===

    // The surface on which the algorithm operates
    HalfedgeMesh *mesh;
    Geometry<Euclidean> *geometry;

    // Internal methods
    void buildOperators();
    VertexData<double> FMMSignedDistance(HalfedgeMesh* mesh, const std::vector<std::pair<VertexPtr, double>>& initialDistances,
                                         const EdgeData<double>& edgeLengths, const HalfedgeData<double>& oppAngles);

    // Small helpers
    GC::SparseMatrix<double> clampedDiagonalInverse(GC::SparseMatrix<double>& A, double thresh);
    VertexData<Eigen::VectorXd> fixed2Dynamic(const VertexData<LabelVec>& input);
    double computeCrossing(double sourceHigh, double sourceLow, double targetHigh, double targetLow);
    double shortestDistanceFromLine(Vector3 target, Vector3 p1, Vector3 p2);
    Vector3 triplePointBarycentrics(Vector3 A, Vector3 B, Vector3 C);
    CornerData<Vector2> parameterize(Geometry<Euclidean>* geom);
    void deleteGeometryAndMesh(Geometry<Euclidean>* &geom);
    void invalidValuePanic(std::string sourceName);


    // Boundary segment helpers
    Vector3 startPoint(const BoundarySegment &b);
    Vector3 endPoint(const BoundarySegment &b);
    double boundaryLength(const BoundarySegment &b);

    // Cached geometric attributes
    void cacheGeometry(void);
    size_t nVert;
    size_t nEdge;
    size_t nFace;
    double surfaceArea = -1;
    double meanEdgeLength = -1;
    double lengthScale = -1;
    double EPS_POSITION;
    size_t INVALID_IND;
    VertexData<size_t> vInd;
    FaceData<size_t> fInd;
    EdgeData<size_t> eInd;
    HalfedgeData<double> cotanWeights;
    EdgeData<double> cotanEdgeWeights;
    EdgeData<double> edgeArea;
    HalfedgeData<double> orientationSign;
    HalfedgeData<double> oppAngles;
    EdgeData<double> edgeLengths;
    DualFaceData<double> vertexArea;
    FaceData<double> faceArea; 
    FaceData<Vector3> faceNormal;
    FaceData<double> faceCurvature; // pointwise, not integrated
    HalfedgeData<Vector3> halfedgeVector;
    VertexData<double> featureSize;

    // Length scale vectors
    GC::DenseVector<double> getScaleVector(bool local=false, int pow=1);
    GC::DenseVector<double> getCutMeshScaleVector(bool local=false, int pow=1);


    // Operators
    GC::SparseMatrix<double> hodge0, hodge1, hodge2, hodge1Inv, hodge0Inv;
    GC::SparseMatrix<double> d0, d0T, d1, d1T;
    GC::SparseMatrix<double> eyeV;

    GC::SparseMatrix<double> zeroFormLaplacian;
    GC::SparseMatrix<double> zeroFormWeakLaplacian;
    GC::SparseMatrix<double> oneFormLaplacian;
    GC::SparseMatrix<double> indicatorSmoother;

    // Used to efficiently zero rows in oneFormLaplacian
    std::vector<std::vector<size_t>> oneFormLaplacianNonzeros;

    // Cache an implicit matrix
    void updateExplicitOperator();
    double cachedLengthCoef = std::numeric_limits<double>::quiet_NaN();
    double cachedBilapCoef = std::numeric_limits<double>::quiet_NaN();
    double cachedVisibilityCoef = std::numeric_limits<double>::quiet_NaN();
    GC::SparseMatrix<double> implicitMat;
    GC::SparseMatrix<double> explicitMat;
    std::unordered_map<double, GC::SparseMatrix<double>> implicitMatCache;

    std::vector<double> interpolateToCutMesh(VertexData<double> vals);
    std::vector<LabelVec> interpolateToCutMesh(VertexData<LabelVec> vals);
};
