// Author: Rohan Sawhney
// Fall 2016

#pragma once

#include "geometry.h"
#include <unordered_map>

typedef pair<double, unsigned int> wvPair;
typedef pair<unsigned int, unsigned int> vvPair;

struct PairHash
{
    static void hashCombine(size_t& seed, const unsigned int& v) {
        hash<unsigned int> hasher;
        seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    
    size_t operator() (const vvPair& p) const {
        size_t hash = 0;
        hashCombine(hash, p.first);
        hashCombine(hash, p.second);
        
        return hash;
    }
};

class MeshCutter {
public:
    // Constructor
    MeshCutter(Geometry<Euclidean> *geometry_);
    
    // Cut
    void cut(const vector<unsigned int>& terminalVertices);
    
protected:
    // Finds shortest path
    double findShortestPath(const unsigned int& src, const unsigned int& target);
    
    // Builds minimal spanning tree
    void buildMinimalSpanningTree(vector<pair<double, vvPair>>& edges,
                                  const vector<unsigned int>& terminalVertices);
    
    // Member variables
    Geometry<Euclidean> *geometry;
    HalfedgeMesh *mesh;
    unordered_map<vvPair, vector<EdgePtr>, PairHash> paths;
    vector<vvPair> spanningTree;
};
