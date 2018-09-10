#include "d_charts.h"

#include <queue>
#include <utility>

DCharts::DCharts(Geometry<Euclidean> *geometry_)
    : geometry(geometry_)
{
    mesh = geometry->getMesh();
    cacheGeometry();
}

void DCharts::cacheGeometry()
{
    geometry->getFaceNormals(normals);
    geometry->getEdgeLengths(edgeLengths);
    geometry->getFaceAreas(faceAreas);

    barycenters = FaceData<Vector3>(mesh);
    for (FacePtr f : mesh->faces())
    {
        barycenters[f] = geometry->barycenter(f);
    }
}

void DCharts::generateCharts()
{

    cout << endl
         << "== Running D-Charts algorithm to generate charts" << endl;

    initializeCharts();

    // Lloyd iterations
    for (size_t iIter = 0; iIter < maxIterations; iIter++)
    {
        cout << "Iteration " << iIter << endl;

        FaceData<int> newLabels = growCharts();

        // Terminate if few enough charts changed
        size_t nChanged = 0;
        for (FacePtr f : mesh->faces())
        {
            if (newLabels[f] != chartLabels[f])
            {
                nChanged++;
            }
        }
        chartLabels = newLabels;

        // Termination criteria
        cout << "Fraction of labels changed: " << ((double)nChanged / mesh->nFaces()) << endl;
        if ((double)nChanged / mesh->nFaces() < terminationFraction)
        {
            cout << "Terminating expansion after " << (iIter + 1) << " iterations due to convergence." << endl;
            break;
        }

        if (iIter + 1 == maxIterations)
        {
            cout << "Terminating expansion after " << maxIterations << " iterations due to iteration limit." << endl;
        }

        // Estimate new cones and centers
        placeChartCenters();
        computeChartCones();
    }
}

void DCharts::initializeCharts()
{

    cout << "Initializing charts" << endl;

    // === Iteratively find fartest points as chart centers
    chartCenters = {mesh->face(0)};
    typedef std::pair<double, FacePtr> FacePair;
    std::priority_queue<FacePair, std::vector<FacePair>, std::greater<FacePair>> faceQueue;
    FacePtr lastFace;
    while (chartCenters.size() < (size_t)numCharts)
    {

        FaceData<double> faceDist(mesh, -1);

        for (FacePtr f : chartCenters)
        {
            faceQueue.push(std::make_pair(0, f));
        }

        while (!faceQueue.empty())
        {

            FacePair curr = faceQueue.top();
            faceQueue.pop();
            double currDist = std::get<0>(curr);
            FacePtr currFace = std::get<1>(curr);

            // Discard stale
            if (faceDist[currFace] != -1)
            {
                continue;
            }

            // Accept
            lastFace = currFace;
            faceDist[currFace] = currDist;

            // Add neighbors
            Vector3 fPos = barycenters[currFace];
            for (FacePtr fNeigh : currFace.adjacentFaces())
            {
                if (faceDist[fNeigh] == -1)
                {
                    Vector3 nPos = barycenters[fNeigh];
                    double nDist = norm(fPos - nPos);
                    faceQueue.push(std::make_pair(currDist + nDist, fNeigh));
                }
            }
        }

        // The last face explored by the search is a new center
        chartCenters.push_back(lastFace);
    }

    // == Assign faces to their closest chart, just for the sake of an initial
    //    computation of the chart cones
    typedef std::pair<double, std::pair<FacePtr, int>> FaceLabelPair;
    std::priority_queue<FaceLabelPair, std::vector<FaceLabelPair>, std::greater<FaceLabelPair>> faceLabelQueue;
    chartLabels = FaceData<int>(mesh, -1);

    // Initialize from seeds
    for (int iC = 0; iC < numCharts; iC++)
    {
        faceLabelQueue.push(std::make_pair(0, std::make_pair(chartCenters[iC], iC)));
    }

    // Grow outwards
    while (!faceLabelQueue.empty())
    {

        FaceLabelPair curr = faceLabelQueue.top();
        faceLabelQueue.pop();
        double currDist = std::get<0>(curr);
        FacePtr currFace = std::get<0>(std::get<1>(curr));
        int currLabel = std::get<1>(std::get<1>(curr));

        // Discard stale
        if (chartLabels[currFace] != -1)
        {
            continue;
        }

        // Accept
        chartLabels[currFace] = currLabel;

        // Add neighbors
        Vector3 fPos = barycenters[currFace];
        for (FacePtr fNeigh : currFace.adjacentFaces())
        {
            if (chartLabels[fNeigh] == -1)
            {
                Vector3 nPos = barycenters[fNeigh];
                double nDist = norm(fPos - nPos);
                faceLabelQueue.push(std::make_pair(currDist + nDist, std::make_pair(fNeigh, currLabel)));
            }
        }
    }

    chartAngles.resize(numCharts);
    chartNormals.resize(numCharts);

    for (int iC = 0; iC < numCharts; iC++)
    {
        chartAngles[iC] = 0.0;
        chartNormals[iC] = normals[chartCenters[iC]];
    }

    // TODO This is different from the initialization described in the paper
    computeChartCones();
}

FaceData<int> DCharts::growCharts()
{
    // Accumulate result here
    FaceData<int> currLabels(mesh, -1);

    // Maintain data about charts
    std::vector<double> chartAreas(numCharts, 0);
    FaceData<double> chartDist(mesh, -1);

    // Compute the cost of adding a face to a chart, except the 1/(chart area) factor
    auto computeBasicCost = [&](FacePtr f, int iChart) {

        double fitError = fittingError(f, iChart);

        double lOuter = 0;
        double lInner = 0;
        double dist = std::numeric_limits<double>::infinity();
        for (HalfedgePtr he : f.adjacentHalfedges())
        {
            if (he.edge().isBoundary())
            {
                lInner += edgeLengths[he.edge()];
            }
            else if (currLabels[he.twin().face()] == iChart)
            {
                lInner += edgeLengths[he.edge()];
                dist = std::min(dist, chartDist[he.twin().face()] + norm(barycenters[f] - barycenters[he.twin().face()]));
            }
            else
            {
                lOuter += edgeLengths[he.edge()];
            }
        }

        double straightError = lOuter / lInner;

        double compactError = PI * dist * dist;

        return std::pow(fitError, alpha) * std::pow(compactError, beta) * std::pow(straightError, gamma);
    };

    // Compute cost, including the area term
    auto computeFullCost = [&](double rawCost, int iChart) {
        return rawCost * std::pow(chartAreas[iChart], -beta);
    };

    // Maintain a separate queue for each chart
    typedef std::pair<double, FacePtr> FaceLabelEntry;
    std::vector<std::priority_queue<FaceLabelEntry, std::vector<FaceLabelEntry>, std::greater<FaceLabelEntry>>> chartQueues(numCharts);

    // Track the current cost for each entry, which will be used to discard stale elements from the queue as costs change
    std::vector<FaceData<double>> chartCosts(numCharts);

    // Initialize with seed triangles
    for (int iC = 0; iC < numCharts; iC++)
    {

        chartCosts[iC] = FaceData<double>(mesh);
        FacePtr seedTri = chartCenters[iC];
        currLabels[seedTri] = iC;
        chartAreas[iC] += faceAreas[seedTri];
        chartDist[seedTri] = 0;

        // Add neighbors to search
        for (FacePtr f : seedTri.adjacentFaces())
        {
            if (currLabels[f] == -1)
            {
                double fitError = fittingError(f, iC);
                if (fitError < fMax)
                {
                    double cost = computeBasicCost(f, iC);
                    chartQueues[iC].push(std::make_pair(cost, f));
                    chartCosts[iC][f] = cost;
                }
            }
        }
    }

    // Grow charts until all triangles belong to a chart
    // TODO This is a simplification of the scheme in the paper, where sufficiently expensive triangles are instead seeded in
    // to new charts
    while (true)
    {

        // Find the smallest entry in any queue
        double minCost = std::numeric_limits<double>::infinity();
        int minChart = -1;
        for (int iC = 0; iC < numCharts; iC++)
        {
            if (!chartQueues[iC].empty())
            {

                const FaceLabelEntry &e = chartQueues[iC].top();
                double rawCost = std::get<0>(e);
                double cost = computeFullCost(rawCost, iC);
                if (cost < minCost)
                {
                    minCost = cost;
                    minChart = iC;
                }
            }
        }

        // Exit if all queues have been drained
        if (minChart == -1)
        {
            break;
        }

        // Process the element from the smallest chart
        FaceLabelEntry e = chartQueues[minChart].top();
        chartQueues[minChart].pop();
        double currCost = std::get<0>(e);
        FacePtr currFace = std::get<1>(e);

        // Discard if stale
        if (currLabels[currFace] != -1 || chartCosts[minChart][currFace] != currCost)
        {
            continue;
        }

        // Accept the entry
        currLabels[currFace] = minChart;
        chartAreas[minChart] += faceAreas[currFace];

        // Recompute distance
        chartDist[currFace] = std::numeric_limits<double>::infinity();
        for (FacePtr fAdj : currFace.adjacentFaces())
        {
            if (currLabels[fAdj] == minChart)
            {
                double dist = chartDist[fAdj] + norm(barycenters[currFace] - barycenters[fAdj]);
                chartDist[currFace] = std::min(chartDist[currFace], dist);
            }
        }

        // Add neighbors to search
        for (FacePtr fNeigh : currFace.adjacentFaces())
        {
            if (currLabels[fNeigh] == -1)
            {
                double fitError = fittingError(fNeigh, minChart);
                if (fitError < fMax)
                {
                    double cost = computeBasicCost(fNeigh, minChart);
                    chartQueues[minChart].push(std::make_pair(cost, fNeigh));
                    chartCosts[minChart][fNeigh] = cost;
                }
            }
        }
    }

    // == Form new charts from unused triangles
    // (the fitting error bound fMax means some triangles may be unplaced)
    for (FacePtr f : mesh->faces())
    {
        if (currLabels[f] != -1)
        {
            continue;
        }

        // Create a new chart
        numCharts++;
        chartCenters.push_back(f);
        chartNormals.push_back(normals[f]);
        chartAngles.push_back(0.0);
        int newChartInd = numCharts-1;

        // Grow a new chart from all connected triangles
        currLabels[f] = newChartInd;
        std::vector<FacePtr> toProcess = {f};
        while(!toProcess.empty()) {

            FacePtr currF = toProcess.back();
            toProcess.pop_back();

            for(FacePtr fNeigh : currF.adjacentFaces()) {
                if(currLabels[fNeigh] == -1) {
                    currLabels[fNeigh] = newChartInd;
                    toProcess.push_back(fNeigh);
                }
            }
        }
    }

    return currLabels;
}

void DCharts::computeChartCones()
{

    for (int iC = 0; iC < numCharts; iC++)
    {

        std::vector<double> chartFaceAreas;
        std::vector<Vector3> chartFaceNormals;
        double totalArea = 0;
        for (FacePtr f : mesh->faces())
        {
            if (chartLabels[f] == iC)
            {
                chartFaceAreas.push_back(faceAreas[f]);
                totalArea += faceAreas[f];
                chartFaceNormals.push_back(normals[f]);
            }
        }

        // Divide out area so we don't have to think about it
        for (size_t i = 0; i < chartFaceAreas.size(); i++)
        {
            chartFaceAreas[i] /= totalArea;
        }

        // Initialize with previous values for chart
        Vector3 normal = chartNormals[iC];
        double cAngle = std::cos(chartAngles[iC]);

        // Compute initial energy
        double oldEnergy = 0;
        for (size_t i = 0; i < chartFaceAreas.size(); i++)
        {
            double error = dot(normal, chartFaceNormals[i]) - cAngle;
            oldEnergy += chartFaceAreas[i] * error * error;
        }

        // === Very simple gradient descent
        double stepSize = 0.2;
        for (int iIter = 0; iIter < 50; iIter++)
        {

            // Compute gradient
            Vector3 gradN{0.0, 0.0, 0.0};
            double gradAngle = 0;
            for (size_t i = 0; i < chartFaceAreas.size(); i++)
            {
                double error = dot(normal, chartFaceNormals[i]) - cAngle;
                gradAngle += -2 * chartFaceAreas[i] * error;
                gradN += 2 * chartFaceAreas[i] * error * chartFaceNormals[i];
            }

            Vector3 newNormal = unit(normal - stepSize * gradN);
            double newCAngle = cAngle - stepSize * gradAngle;

            // Check the energy of the step
            double energy = 0;
            for (size_t i = 0; i < chartFaceAreas.size(); i++)
            {
                double error = dot(newNormal, chartFaceNormals[i]) - newCAngle;
                energy += chartFaceAreas[i] * error * error;
            }

            // Check convergence and tweak step size
            if (std::abs(angle(normal, newNormal)) < 1e-4 && std::abs(cAngle - newCAngle) < 1e-4)
            {
                break;
            }
            if (energy > oldEnergy)
            {
                stepSize *= 0.5;
            }
            else
            {
                // Step and project
                normal = newNormal;
                cAngle = newCAngle;
                oldEnergy = energy;
            }
        }

        chartNormals[iC] = normal;
        chartAngles[iC] = std::acos(cAngle);
    }
}

void DCharts::placeChartCenters()
{

    // TODO I'm not acutally sure what the proper way to do this is...
    // This method diverges from the paper, and simply uses the face which is most distant from
    // any part of the boundary

    for (int iC = 0; iC < numCharts; iC++)
    {

        typedef std::pair<double, FacePtr> FaceDistEntry;
        std::priority_queue<FaceDistEntry, std::vector<FaceDistEntry>, std::greater<FaceDistEntry>> queue;
        FaceData<double> faceDist(mesh, -1);

        // Initialize search with boundary faces
        FacePtr lastFace;
        for (FacePtr f : mesh->faces())
        {
            if (chartLabels[f] == iC)
            {
                bool isBoundary = false;
                for (FacePtr fAdj : f.adjacentFaces())
                {
                    if (chartLabels[fAdj] != iC)
                    {
                        isBoundary = true;
                    }
                }

                if (isBoundary)
                {
                    faceDist[f] = 0.0;
                    lastFace = f;

                    // Add neighbors
                    for (FacePtr fAdj : f.adjacentFaces())
                    {
                        if (chartLabels[fAdj] == iC)
                        {
                            double dist = norm(barycenters[f] - barycenters[fAdj]);
                            queue.push(std::make_pair(dist, fAdj));
                        }
                    }
                }
            }
        }

        // Search inward from boundary
        while (!queue.empty())
        {

            FaceDistEntry e = queue.top();
            queue.pop();
            double currDist = std::get<0>(e);
            FacePtr currFace = std::get<1>(e);

            // Discard stale
            if (faceDist[currFace] != -1)
            {
                continue;
            }

            // Accept
            faceDist[currFace] = currDist;
            lastFace = currFace;

            // Add neighbors
            for (FacePtr fAdj : currFace.adjacentFaces())
            {
                if (chartLabels[fAdj] == iC)
                {
                    double dist = currDist + norm(barycenters[currFace] - barycenters[fAdj]);
                    queue.push(std::make_pair(dist, fAdj));
                }
            }
        }

        // The seed is the last face we found; the most distant from the boundary
        chartCenters[iC] = lastFace;
    }
}

double DCharts::fittingError(FacePtr f, int iChart)
{
    double e = dot(chartNormals[iChart], normals[f]) - std::cos(chartAngles[iChart]);
    return e * e;
}