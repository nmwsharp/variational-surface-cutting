#pragma once

#include "geometry.h"

class DCharts
{

  public:
    DCharts(Geometry<Euclidean> *geometry_);

    // Parameters
    int numCharts = 10;
    double alpha = 1.0; // as in paper
    double beta = 0.7;  // as in paper
    double gamma = 0.5; // as in paper
    // double alpha = 1.0; // as in paper
    // double beta = 0.0;  // as in paper
    // double gamma = 0.0; // as in paper
    size_t maxIterations = 500;
    double terminationFraction = 0.1;
    double fMax = 0.2;

    // State
    FaceData<int> chartLabels; // solution is here
    std::vector<Vector3> chartNormals;
    std::vector<double> chartAngles;
    std::vector<FacePtr> chartCenters;

    // Runs all of the steps below until termiation, result is stored in chart labels
    void generateCharts();

    void estimateNumCharts();
    void initializeCharts();
    FaceData<int> growCharts();
    void computeChartCones();
    void placeChartCenters();

  private:
    HalfedgeMesh *mesh;
    Geometry<Euclidean> *geometry;

    void cacheGeometry();
    FaceData<Vector3> normals;
    FaceData<Vector3> barycenters;
    FaceData<double> faceAreas;
    EdgeData<double> edgeLengths;

    // Small helpers
    double fittingError(FacePtr f, int iChart);
};