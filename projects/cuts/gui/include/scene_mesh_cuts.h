#pragma once

#include "scene_mesh.h"
#include "scene_vectors.h"
#include "scene_lines.h"
#include "detect_symmetry.h"

class SceneMeshCuts : public SceneMesh {

    public:

        SceneMeshCuts(Viewer& parent_, Geometry<Euclidean>* geometryToDraw);
        ~SceneMeshCuts();

        virtual void draw() override;

        // === On the shape
        void showRegionLabels(const FaceData<int>& regionLabels);
        void showRegionLabels(const VertexData<Eigen::VectorXd>& vals);
        void showRegionChecker(FaceData<int> regionLabels, CornerData<Vector2> params);
        void showRegionChecker(CornerData<Vector2> params);
        void showScalar(const VertexData<double>& vals, bool symmetric = false);
        void showScalar(const FaceData<double>& vals, bool symmetric = false);
        void showMultiScalar(const FaceData<Eigen::VectorXd>& vals);
        void showMultiScalar(const VertexData<Eigen::VectorXd>& vals);
        void showSymmetry(const SymmetryResult& s);
        void showNothing();

        void showVector(const FaceData<Vector3>& vals);
        void showOneForm(const EdgeData<double>& oneForm, bool sparse);
        void showEdgePerpVector(const EdgeData<double>& vals);
        void showEdgeVector(const EdgeData<double>& vals);
        void showNothingVector();


        // === Other 
        void setLines(const std::vector<std::array<Vector3,2>>& lines);
        void toggleLines(bool b);
        void setLinesAlternate(const std::vector<std::array<Vector3,2>>& lines);
        void toggleLinesAlternate(bool b);

        void setCurrentMesh(Geometry<Euclidean>* geometry);

    private:

        void computeSymmetricBounds(VertexData<double> data, double& min, double &max, double cap = std::numeric_limits<double>::infinity());
        void computeSymmetricBounds(FaceData<double> data, double& min, double &max, double cap = std::numeric_limits<double>::infinity());

        double shapeLengthScale = -1;
        SceneVectors vectorArtist;
        SceneLines lineArtist;
        SceneLines lineArtistAlternate;

        // Colors to use
        // std::vector<Vector3> palette = {
        //                             Vector3{244,67,54},   // red
        //                             Vector3{63,81,181},   // indigo
        //                             Vector3{0,150,136},   // teal
        //                             Vector3{255,152,0},   // orange
        //                             Vector3{156,39,176},  // purple
        //                             };

        Vector3 getPaletteColor(size_t iRegion);
};
