// Author: Rohan Sawhney
// Fall 2016

#include "mesh_cutter.h"
#include <queue>
#include <functional>

#include "disjoint_sets.h"

MeshCutter::MeshCutter(Geometry<Euclidean> *geometry_):
geometry(geometry_),
mesh(geometry->getMesh())
{
    
}

double MeshCutter::findShortestPath(const unsigned int& src, const unsigned int& target)
{
    // Initialize
    VertexData<size_t> vIndices = mesh->getVertexIndices();
    HalfedgeData<size_t> hIndices = mesh->getHalfedgeIndices();
    
    vvPair stPair(src, target);
    bool boundaryVerts = mesh->vertex(src).isBoundary() && mesh->vertex(target).isBoundary();
    unsigned int vCount = (int)mesh->nVertices();
    unsigned int hReal = (int)mesh->nHalfedges();
    unsigned int hCount = hReal + mesh->nImaginaryHalfedges();
    
    priority_queue<wvPair, vector<wvPair>, greater<wvPair>> pq;
    vector<double> dist(vCount, INFINITY);
    vector<unsigned int> prev(vCount, hCount);
    
    // Push root
    pq.push(make_pair(0.0, src));
    dist[src] = 0.0;
    
    // Loop till queue becomes empty
    while (!pq.empty()) {
        // Get vertex on top of queue
        unsigned int uIdx = pq.top().second;
        pq.pop();
        
        VertexPtr u = mesh->vertex(uIdx);
        
        // Terminate search if target has been found
        if (uIdx == target) break;
        
        for (HalfedgePtr h: u.incomingHalfedges()) {
            // Get neighboring vertex index
            VertexPtr v = h.vertex();
            unsigned int vIdx = vIndices[v];
            
            // Compute edge weight
            EdgePtr e = h.edge();
            double weight = boundaryVerts && e.isBoundary() ?
                            numeric_limits<double>::min() : geometry->length(e);
            
            // Update distance
            if (dist[vIdx] > dist[uIdx] + weight) {
                dist[vIdx] = dist[uIdx] + weight;
                prev[vIdx] = hIndices[h.twin()];
                pq.push(make_pair(dist[vIdx], vIdx));
            }
        }
    }
    
    // Construct path
    double weight = 0.0;
    unsigned int uIdx = target;
    while (prev[uIdx] != hCount) {
        unsigned int hIdx = prev[uIdx];
        HalfedgePtr h = hIdx < hReal ? mesh->halfedge(hIdx) : mesh->imaginaryHalfedge(hIdx - hReal);
        EdgePtr e = h.edge();
        
        paths[stPair].push_back(e);
        weight += boundaryVerts && e.isBoundary() ?
                  numeric_limits<double>::min() : geometry->length(e);
        uIdx = vIndices[h.vertex()];
    }
    
    return weight;
}

void MeshCutter::buildMinimalSpanningTree(vector<pair<double, vvPair>>& edges,
                                          const vector<unsigned int>& terminalVertices)
{
    // Sort
    sort(edges.begin(), edges.end());
    
    // Create disjoint sets
    DisjointSets ds(edges.size());
    
    // Iterate through all sorted edges
    for (vector<pair<double, vvPair>>::const_iterator it = edges.begin(); it != edges.end(); it++) {
        unsigned int uIdx = it->second.first;
        unsigned int vIdx = it->second.second;
        
        unsigned int setU = ds.find(uIdx);
        unsigned int setV = ds.find(vIdx);
        
        // If edge is not creating a cycle then add it to the MST
        if (setU != setV) {
            unsigned int v1 = terminalVertices[uIdx];
            unsigned int v2 = terminalVertices[vIdx];
            if (v1 > v2) swap(v1, v2);
            spanningTree.push_back(make_pair(v1, v2));
            
            // Merge sets
            ds.merge(setU, setV);
        }
    }
}

void MeshCutter::cut(const vector<unsigned int>& terminalVertices)
{
    // Unmark any previouly existing cut edges
    for (EdgePtr e: mesh->edges()) e.markCut(false);
    
    // Find shortest paths
    vector<pair<double, vvPair>> edges;
    for (unsigned int i = 0; i < (unsigned int)terminalVertices.size()-1; i++) {
        for (unsigned int j = i+1; j < (unsigned int)terminalVertices.size(); j++) {
            unsigned int v1 = terminalVertices[i];
            unsigned int v2 = terminalVertices[j];
            if (v1 > v2) swap(v1, v2);
            
            double weight = findShortestPath(v1, v2);
            edges.push_back(make_pair(weight, make_pair(i, j)));
        }
    }
    
    // Build minimal spanning tree
    buildMinimalSpanningTree(edges, terminalVertices);
    
    // Mark cut edges
    for (int i = 0; i < (int)spanningTree.size(); i++) {
        vvPair& p = spanningTree[i];
        for (int j = 0; j < (int)paths[p].size(); j++) {
            EdgePtr e = paths[p][j];
            e.markCut(true);
        }
    }
}
