// Author: Rohan Sawhney
// Fall 2016

#include "flatten_utils.h"

bool UnconstraintedMinimization::solve(GC::DenseMatrix<double>& x, Handle& handle) {
    int n = x.length();
    int k = 0;
    double f = 0.0;
    handle.computeEnergy(f, x);
    
    while (true) {
        // Compute update direction
        GC::DenseMatrix<double> g(n);
        handle.computeGradient(g, x);
        
        GC::SparseMatrix<double> H(n, n);
        handle.computeHessian(H, x);
        
        GC::DenseMatrix<double> p(n);
        solvePositiveDefinite(H, p, g);
        
        // Compute step size
        double t = 1.0;
        double fp = f;
        handle.computeEnergy(f, x - t*p);
        while (f > fp - 0.5*t*dot(g, p)) {
            t = 0.9*t;
            handle.computeEnergy(f, x - t*p);
        }
        
        // Terminate if f is not finite
        if (!isfinite(f)) return false;
        
        // Update
        x -= t*p;
        k++;
        
        // Check termination condition
        if (g.norm() < 1e-8 || fabs(f - fp) < 1e-8 || k > 1000) break;
    }
    
    return true;
}

Layout::AngleLengthLayout::AngleLengthLayout(HalfedgeMesh *mesh_,
                                             const CornerData<size_t>& wIndices_,
                                             const CornerData<double>& angles_,
                                             const EdgeData<double>& lengths_):
mesh(mesh_),
fIndices(mesh->getFaceIndices()),
wIndices(wIndices_),
angles(angles_),
lengths(lengths_)
{
    
}

void Layout::AngleLengthLayout::performFaceLayout(HalfedgePtr h, const Vector2& dir, vector<Vector2>& x,
                                                  unordered_map<int, bool>& visited, stack<EdgePtr>& stack)
{
    int fIdx = fIndices[h.face()];
    if (visited.find(fIdx) == visited.end()) {
        HalfedgePtr next = h.next();
        HalfedgePtr prev = h.prev();
        
        // Compute new direction from angle
        double angle = angles[next.corner()];
        Vector2 newDir = {cos(angle)*dir[0] - sin(angle)*dir[1],
                          sin(angle)*dir[0] + cos(angle)*dir[1]};
        
        // Compute position
        x[wIndices[h.corner()]] = x[wIndices[next.corner()]] + newDir*lengths[prev.edge()];
        
        // Mark face as visited
        visited[fIdx] = true;
        
        // Push edges onto stack
        stack.push(next.edge());
        stack.push(prev.edge());
    }
}

void Layout::AngleLengthLayout::performLayout(vector<Vector2>& x)
{
    // Push any interior edge on the stack
    stack<EdgePtr> stack;
    for (EdgePtr e: mesh->edges()) {
        if (!e.isBoundary() && !e.isCut()) {
            stack.push(e);
            break;
        }
    }
    
    // Place this edge on the x axis
    EdgePtr e = stack.top();
    HalfedgePtr h = e.halfedge();
    unsigned int w1 = wIndices[h.next().corner()];
    unsigned int w2 = wIndices[h.prev().corner()];
    
    x[w1] = Vector2{0, 0};
    x[w2] = Vector2{lengths[e], 0};
    
    // Perform layout
    unordered_map<int, bool> visited;
    while (!stack.empty()) {
        e = stack.top();
        stack.pop();
        
        // Only perform layout from interior edges
        if (!e.isBoundary() && !e.isCut()) {
            // Compute halfedge direction
            h = e.halfedge();
            w1 = wIndices[h.next().corner()];
            w2 = wIndices[h.prev().corner()];
            Vector2 dir = x[w2] - x[w1];
            dir.normalize();
            
            // Perform face layout for halfedge and its twin
            performFaceLayout(h, dir, x, visited, stack);
            performFaceLayout(h.twin(), -dir, x, visited, stack);
        }
    }
}
