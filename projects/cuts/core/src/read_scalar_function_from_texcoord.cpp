#include "read_scalar_function_from_texcoord.h"

#include <vector>
#include <sstream>
#include <fstream>

namespace {

class Index {
public:
    Index() {}
    
    Index(int v, int vt, int vn) : position(v), uv(vt), normal(vn) {}
    
    bool operator<(const Index& i) const {
        if (position < i.position) return true;
        if (position > i.position) return false;
        if (uv < i.uv) return true;
        if (uv > i.uv) return false;
        if (normal < i.normal) return true;
        if (normal > i.normal) return false;
        
        return false;
    }
    
    int position;
    int uv;
    int normal;
};

Index parseFaceIndex(const std::string& token)
{
    std::stringstream in(token);
    std::string indexString;
    int indices[3] = {-1, -1, -1};
    
    int i = 0;
    while (std::getline(in, indexString, '/')) {
        if (indexString != "\\") {
            std::stringstream ss(indexString);
            ss >> indices[i++];
        }
    }
    
    // decrement since indices in OBJ files are 1-based
    return Index(indices[0]-1,
                 indices[1]-1,
                 indices[2]-1);
}

// Read a .obj file containing a polygon mesh
std::vector<double> parseOBJForScalarAtVerts(std::string filename)
{
    std::cout << "Reading mesh from file: "<< filename << std::endl;

    std::vector<double> texCoordVals;
    std::vector<double> vertexVals;
    
	// Open the file
	std::ifstream in(filename);
	if (!in) throw std::invalid_argument("Could not open mesh file " + filename);
    
    // parse obj format
    std::string line;
    while (getline(in, line)) {
        std::stringstream ss(line);
        std::string token;
        
        ss >> token;
        
        if (token == "v") {
            
        } else if (token == "vt") {
            double u, v;
            ss >> u >> v;
            
            texCoordVals.push_back(u);
            
        } else if (token == "vn") {
            // Do nothing
            
        } else if (token == "f") {
            while (ss >> token) {
                Index index = parseFaceIndex(token);
                if (index.position < 0) {
                    getline(in, line);
                    size_t i = line.find_first_not_of("\t\n\v\f\r ");
                    index = parseFaceIndex(line.substr(i));
                }

                // Copy the function to vertices
                int iVertex = index.position;
                int iTex = index.uv;

                // Make sure the list of vertices is long enough
                if(iVertex >= (int)vertexVals.size()) {
                    if(vertexVals.capacity() < texCoordVals.size()) {
                        vertexVals.reserve(texCoordVals.size()); // reserve extra to avoid O(N^2) expansion
                    }
                    vertexVals.resize(iVertex+1);
                }

                vertexVals[iVertex] = texCoordVals[iTex];
            
            }
            
        }
    }

    return vertexVals;
}
} // Anonymous namespace

// Assumptions: objFile should correspond to mesh, and vertex order
// must match
// Scalar function should be stored at 'u' coordinate of texcoords,
// and must take a single value around each vetex
VertexData<double> readScalarFromTex(HalfedgeMesh* mesh, std::string objFile) {

    // Parse the values out of the texture coordinates
    std::vector<double> vertVals = parseOBJForScalarAtVerts(objFile);

    if(vertVals.size() != mesh->nVertices()) {
        throw std::runtime_error("Number of vertices in mesh (" + std::to_string(mesh->nVertices()) + ") does not match number from obj (" + std::to_string(vertVals.size()) + ")");
    }

    // Transfer to vertex data
    VertexData<double> result(mesh);
    for(size_t i = 0; i < mesh->nVertices(); i++) {
        // result[mesh->vertex(i)] = vertVals[i] * (1.0 / .003908); // undo stupid scaling on export
        result[mesh->vertex(i)] = vertVals[i];
    }

    return result;
}