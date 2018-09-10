#include<scene_mesh.h>

#include <iostream>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <limits>
 
#include <shaders.h>
#include <shaders/shiny_shaders.h>
#include <shaders/phong_shaders.h>
#include <shaders/wire_shaders.h>
#include <shaders/light_shaders.h>
#include <shaders/shadow_shaders.h>
#include <shaders/dual_shaders.h>
#include <shaders/checker_shaders.h>
#include <shaders/simple_shaders.h>

using std::cout; using std::endl;

SceneMesh::SceneMesh(Viewer &parent_, Geometry<Euclidean>* geometryToDraw):
        SceneObject(parent_),
        faceGLProgram(new GLProgram(&SHINY_VERT_SHADER, &SHINY_FRAG_SHADER, DrawMode::Triangles)),
        //faceGLProgram(&PHONG_VERT_SHADER, &PHONG_FRAG_SHADER, DrawMode::Triangles),
        wireGLProgram(new GLProgram(&WIRE_VERT_SHADER, &WIRE_FRAG_SHADER, DrawMode::Triangles)),
        dualGLProgram(new GLProgram(&PASSTHRU_DUAL_VERT_SHADER, &DUAL_GEOM_SHADER, &DUAL_FRAG_SHADER, DrawMode::LinesAdjacency)),
        lightGLProgram(new GLProgram(&LIGHT_VERT_SHADER, &LIGHT_FRAG_SHADER, DrawMode::Triangles)),
        //shadowGLProgram(&PASSTHRU_SHADOW_VERT_SHADER, &SHADOW_GEOM_SHADER, &SHADOW_FRAG_SHADER, DrawMode::TrianglesAdjacency),
        shadowGLProgram(nullptr),
        pickGLProgram(new GLProgram(&SIMPLE_VERT_SHADER, &SIMPLE_FRAG_SHADER, DrawMode::Triangles)),
        mesh(geometryToDraw->getMesh()),
        geometry(geometryToDraw)
{
   // Initialize program for shadow volumes
   if(mesh->isSimplicial()) {
      // Triangle adjacency is easy to infer from a pure triangle
      // mesh, so we use the standard "silhouettes at edges" approach.
      shadowGLProgram = new GLProgram( &SHADOW_VERT_SHADER,
                                       &SHADOW_GEOM_SHADER_TRIADJ,
                                       &SHADOW_FRAG_SHADER,
                                       DrawMode::TrianglesAdjacency );
   } else {
      // For general polygonal meshes, we locally tessellate each
      // polygon, making adjacency more difficult to infer.  Hence we
      // use the silhouette of the interpolated vertex normals, which
      // can be computed from a single triangle (but tends to produce
      // more artifacts...)
      shadowGLProgram = new GLProgram( &SHADOW_VERT_SHADER,
                                       &SHADOW_GEOM_SHADER_TRI,
                                       &SHADOW_FRAG_SHADER,
                                       DrawMode::Triangles );
   }

    // Compute a center and scale for the data, then adjust the viewer's camera accordingly
    dataCenter = geometry->center();
    dataScale = geometry->lengthScale();

    // Fill GL buffers for drawing the geometry
    fillMeshBuffers();
    //fillMeshBuffersDual();
    fillMeshBuffersShadow();
    fillMeshBuffersPick();
}

SceneMesh::~SceneMesh() {
    if (faceGLProgram != nullptr) delete faceGLProgram;
    if (wireGLProgram!= nullptr) delete wireGLProgram;
    if(dualGLProgram != nullptr) delete dualGLProgram;
    if(lightGLProgram != nullptr) delete lightGLProgram;
    if(shadowGLProgram != nullptr) delete shadowGLProgram;
    if(pickGLProgram != nullptr) delete pickGLProgram;
}

void SceneMesh::draw() {

    // Set uniforms
    glm::mat4 viewMat = parent.camera.getViewMatrix();
    faceGLProgram->setUniform("u_viewMatrix", glm::value_ptr(viewMat));
    wireGLProgram->setUniform("u_viewMatrix", glm::value_ptr(viewMat));
    dualGLProgram->setUniform("u_viewMatrix", glm::value_ptr(viewMat));

    glm::mat4 projMat = parent.camera.getPerspectiveMatrix();
    faceGLProgram->setUniform("u_projMatrix", glm::value_ptr(projMat));
    wireGLProgram->setUniform("u_projMatrix", glm::value_ptr(projMat));
    dualGLProgram->setUniform("u_projMatrix", glm::value_ptr(projMat));

    Vector3 eyePos = parent.camera.getCameraWorldPosition();
    faceGLProgram->setUniform("u_eye", eyePos);

    // The light is over the right shoulder of the camera
    //Vector3 zAxis{0.0, 0.0, 1.0};
    //Vector3 lightPos = 2*(eyePos.rotate_around(zAxis, 3.14159 / 8) + Vector3{0.0, 0.0, norm(eyePos)/5.0});
    Vector3 lightPos = parent.camera.getLightWorldPosition();
    faceGLProgram->setUniform("u_light", lightPos);

    wireGLProgram->setUniform("u_wireframeWeight", parent.wireframeOpacity);
    
    Vector3 wireColor{ 0., 0., 0. };
    wireGLProgram->setUniform("u_wireColor", wireColor);
    
    glEnable( GL_DEPTH_TEST );
    glDisable( GL_BLEND );
    faceGLProgram->draw();

    if(parent.showEdges) {
       glEnable(GL_BLEND);
       glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
       glDepthMask(GL_FALSE);
       glDepthFunc(GL_LEQUAL);
       wireGLProgram->draw();
       glDepthMask(GL_TRUE);
       glDisable(GL_BLEND);
    }

    // XXX no longer works for general polygon meshes
    // if(parent.showDual) {
    //    //glEnable(GL_LINE_SMOOTH);
    //    //glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
    //    //glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    //    glDepthMask(GL_FALSE);
    //    dualGLProgram->draw();
    //    glDepthMask(GL_TRUE);
    //    //glDisable(GL_LINE_SMOOTH);
    //    //glDisable(GL_BLEND);
    // }
}

void SceneMesh::drawLight() {
   glm::mat4 viewMat = parent.camera.getViewMatrix();
   lightGLProgram->setUniform("u_viewMatrix", glm::value_ptr(viewMat));

   glm::mat4 projMat = parent.camera.getPerspectiveMatrix();
   lightGLProgram->setUniform("u_projMatrix", glm::value_ptr(projMat));

   Vector3 eyePos = parent.camera.getCameraWorldPosition();
   lightGLProgram->setUniform("u_eye", eyePos);

   Vector3 lightPos = parent.camera.getLightWorldPosition();
   lightGLProgram->setUniform("u_light", lightPos);

   lightGLProgram->draw();
}

void SceneMesh::drawShadow() {
   glm::mat4 viewMat = parent.camera.getViewMatrix();
   shadowGLProgram->setUniform("u_viewMatrix", glm::value_ptr(viewMat));

   glm::mat4 projMat = parent.camera.getPerspectiveMatrix();
   shadowGLProgram->setUniform("u_projMatrix", glm::value_ptr(projMat));

   Vector3 lightPos = parent.camera.getLightWorldPosition();
   shadowGLProgram->setUniform("u_light", lightPos);
   shadowGLProgram->draw();
}

void SceneMesh::drawPick()
{
    glm::mat4 viewMat = parent.camera.getViewMatrix();
    pickGLProgram->setUniform("u_viewMatrix", glm::value_ptr(viewMat));
    
    glm::mat4 projMat = parent.camera.getPerspectiveMatrix();
    pickGLProgram->setUniform("u_projMatrix", glm::value_ptr(projMat));
    
    pickGLProgram->draw();
}
    
void SceneMesh::resetCameraForData() {
    parent.camera.dist = dataScale * 0.5;
    parent.camera.dataLengthScale = dataScale;
    parent.camera.lookAtPoint = dataCenter;
}

void SceneMesh::updateGlBuffers(std::vector<Vector3> positions, std::vector<Vector3> normals,
                                std::vector<Vector3> vertexNormals, std::vector<Vector3> barycentric,
                                bool update)
{
    faceGLProgram->setAttribute("a_position", positions, update);
    faceGLProgram->setAttribute("a_normal", normals, update);
    
    wireGLProgram->setAttribute("a_position", positions, update);
    wireGLProgram->setAttribute("a_baryCoord", barycentric, update);
    
    lightGLProgram->setAttribute("a_position", positions, update);
    lightGLProgram->setAttribute("a_normal", normals, update);
    
    if(!mesh->isSimplicial())
    {
        shadowGLProgram->setAttribute("a_position", positions, update);
        shadowGLProgram->setAttribute("a_normal", vertexNormals, update);
    }
}

void SceneMesh::fillMeshBuffers(bool update) {

    bufferMode = FACE_BASED_BUFFERS;
    std::vector<Vector3> positions;
    std::vector<Vector3> faceNormals;
    std::vector<Vector3> vertexNormals;
    std::vector<Vector3> barycentric;
 
    if(mesh->isSimplicial()) { // Pure triangle mesh

       unsigned int i = 0;
       positions.resize(3*mesh->nFaces());
       faceNormals.resize(3*mesh->nFaces());
       vertexNormals.resize(3*mesh->nFaces());
       barycentric.resize(3*mesh->nFaces());
       for(FacePtr f : mesh->faces()) {
          HalfedgePtr he = f.halfedge();
          for(unsigned int j = 0; j < 3; j++) {
             VertexPtr v = he.vertex();
             positions[3*i + j] = geometry->position(v);
             faceNormals[3*i + j] = geometry->normal(f);
             vertexNormals[3*i + j] = geometry->normal(v);
             he = he.next();
          }
          barycentric[3*i+0] = Vector3{1.,0.,0.};
          barycentric[3*i+1] = Vector3{0.,1.,0.};
          barycentric[3*i+2] = Vector3{0.,0.,1.};
          i++;
       }

    } else { // General polygon mesh

       for(FacePtr f : mesh->faces()) {
          if( f.degree() == 3 ) { // triangle
             HalfedgePtr he = f.halfedge();
             for(unsigned int j = 0; j < 3; j++) {
                VertexPtr v = he.vertex();
                positions.push_back( geometry->position(v) );
                faceNormals.push_back( geometry->normal(f) );
                vertexNormals.push_back( geometry->normal(v) );
                he = he.next();
             }
             barycentric.push_back( Vector3{1.,0.,0.} );
             barycentric.push_back( Vector3{0.,1.,0.} );
             barycentric.push_back( Vector3{0.,0.,1.} );
          } else { // polygon

             std::vector<Triangle> triangles = f.triangulation();
             size_t k = triangles.size();
             for(size_t j = 0; j < k; j++ ) {
                for(size_t i=0; i < 3; i++) {
                   VertexPtr v = triangles[j][i];
                   positions.push_back( geometry->position(v) );
                   faceNormals.push_back( geometry->normal(f) );
                   vertexNormals.push_back( geometry->normal(v) );
                }

                // Hide edges that are not part of the original
                // tessellation by setting the opposite barycentric
                // coordinates to "infinity" (these values are
                // used by the shader to compute the edge function).
                const double X = std::numeric_limits<double>::infinity();
                if( j==0 || j==k-1 )
                {
                   barycentric.push_back(Vector3{X,0,0});
                   barycentric.push_back(Vector3{0,1,0});
                   barycentric.push_back(Vector3{0,0,1});
                }
                else
                {
                   barycentric.push_back(Vector3{X,0,0});
                   barycentric.push_back(Vector3{0,X,0});
                   barycentric.push_back(Vector3{0,0,1});
                }
             }
          }
       }
    }

    std::vector<Vector3>& normals = parent.flatShaded ? faceNormals : vertexNormals;
    updateGlBuffers(positions, normals, vertexNormals, barycentric, update);
}

void SceneMesh::fillVertexBasedMeshBuffers(bool update)
{
    bufferMode = VERTEX_BASED_BUFFERS;
    std::vector<Vector3> positions;
    std::vector<Vector3> faceNormals;
    std::vector<Vector3> barycentric;
    
    if (mesh->isSimplicial()) { // Pure triangle mesh
        unsigned int i = 0;
        positions.resize(18*mesh->nFaces());
        faceNormals.resize(18*mesh->nFaces());
        barycentric.resize(18*mesh->nFaces());
        for (FacePtr f : mesh->faces()) {
            Vector3 b = geometry->barycenter(f);
            unsigned int j = 0;
            for (HalfedgePtr h: f.adjacentHalfedges()) {
                Vector3 m = geometry->midpoint(h.edge());
                
                Vector3 v = geometry->position(h.vertex());
                positions[18*i+6*j+0] = v;
                positions[18*i+6*j+1] = m;
                positions[18*i+6*j+2] = b;
                barycentric[18*i+6*j+0] = Vector3{1.,0.,0.};
                barycentric[18*i+6*j+1] = Vector3{0.,1.,0.};
                barycentric[18*i+6*j+2] = Vector3{0.,0.,1.};
                Vector3 n = cross(m-v, b-v); n.normalize();
                for (unsigned int k = 0; k < 3; k++) faceNormals[18*i+6*j+k] = n;
                
                Vector3 vn = geometry->position(h.next().vertex());
                positions[18*i+6*j+3] = b;
                positions[18*i+6*j+4] = m;
                positions[18*i+6*j+5] = vn;
                barycentric[18*i+6*j+3] = Vector3{1.,0.,0.};
                barycentric[18*i+6*j+4] = Vector3{0.,1.,0.};
                barycentric[18*i+6*j+5] = Vector3{0.,0.,1.};
                Vector3 nn = cross(m-b, vn-b); nn.normalize();
                for (unsigned int k = 0; k < 3; k++) faceNormals[18*i+6*j+3+k] = nn;
                
                j++;
            }
            i++;
        }
    
    } else { // General polygon mesh
        // TODO: The flatten tool does not work with general polygon meshes
        // so I did not implement this part.
    }
    
    updateGlBuffers(positions, faceNormals, faceNormals, barycentric, update);
}

void SceneMesh::fillEdgeBasedMeshBuffers(bool update)
{
    bufferMode = EDGE_BASED_BUFFERS;
    std::vector<Vector3> positions;
    std::vector<Vector3> faceNormals;
    std::vector<Vector3> barycentric;
    
    if (mesh->isSimplicial()) { // Pure triangle mesh
        unsigned int i = 0;
        positions.resize(9*mesh->nFaces());
        faceNormals.resize(9*mesh->nFaces());
        barycentric.resize(9*mesh->nFaces());
        for (FacePtr f : mesh->faces()) {
            Vector3 b = geometry->barycenter(f);
            unsigned int j = 0;
            for (HalfedgePtr h: f.adjacentHalfedges()) {
                Vector3 v = geometry->position(h.vertex());
                Vector3 vn = geometry->position(h.next().vertex());
                
                positions[9*i+3*j+0] = b;
                positions[9*i+3*j+1] = v;
                positions[9*i+3*j+2] = vn;
                barycentric[9*i+3*j+0] = Vector3{1.,0.,0.};
                barycentric[9*i+3*j+1] = Vector3{0.,1.,0.};
                barycentric[9*i+3*j+2] = Vector3{0.,0.,1.};
                Vector3 n = cross(v-b, vn-b); n.normalize();
                for (unsigned int k = 0; k < 3; k++) faceNormals[9*i+3*j+k] = n;
                
                j++;
            }
            i++;
        }
        
    } else { // General polygon mesh
        // TODO: The flatten tool does not work with general polygon meshes
        // so I did not implement this part.
    }
    
    updateGlBuffers(positions, faceNormals, faceNormals, barycentric, update);
}

void SceneMesh::fillMeshBuffersShadow(bool update) {

   if(mesh->isSimplicial()) {
      // Traditional edge-based method (using triangles_adjacency)
      std::vector<Vector3> positions(6*mesh->nFaces());

      // Walk the faces building array of positions
      unsigned int i = 0;
      for(FacePtr f : mesh->faces()) {
         HalfedgePtr he = f.halfedge();
         for(unsigned int j = 0; j < 3; j++) {
            VertexPtr v = he.twin().vertex();
            VertexPtr w = he.next().twin().next().twin().vertex();
            positions[6*i + 2*j + 0] = geometry->position(v);
            if(he.next().twin().isReal()) {
               positions[6*i + 2*j + 1] = geometry->position(w);
            } else {
               positions[6*i + 2*j + 1] = geometry->position(v) + 0.5 * geometry->vector(he.next()) - 0.1 * geometry->vector(he);
            }
            he = he.next();
         }
         i++;
      }

      // Store the values in GL buffers
      shadowGLProgram->setAttribute("a_position", positions, update);
      shadowGLProgram->setAttribute("a_normal", positions, update); // not used, so we can just pass a copy of the positions
   }
}

void SceneMesh::fillMeshBuffersDual(bool update)
{
   unsigned int i;
   std::vector<Vector3> positions(4*mesh->nFaces());

   i = 0;
   for(DualVertexPtr v : mesh->dual().vertices()) {
      positions[4*i+0] = geometry->position(v);
      i++;
   }

   i=0;
   for(FacePtr f : mesh->faces()) {
      unsigned int j = 1;
      for(HalfedgePtr h : f.adjacentHalfedges()) {
         Vector3 a = geometry->position(h.vertex());
         Vector3 b = geometry->position(h.twin().vertex());
         positions[4*i+j] = (a+b)/2.;
         j++;
      }
      i++;
   }

   dualGLProgram->setAttribute("a_position", positions, update);
}

Vector3 elementColor(unsigned int idx)
{
    return Vector3{((idx & 0x000000FF) >>  0)/255.0,
                   ((idx & 0x0000FF00) >>  8)/255.0,
                   ((idx & 0x00FF0000) >> 16)/255.0};
}

void SceneMesh::fillMeshBuffersPick(bool update)
{
    VertexData<size_t> vIndices = mesh->getVertexIndices();
    EdgeData<size_t> eIndices = mesh->getEdgeIndices();
    FaceData<size_t> fIndices = mesh->getFaceIndices();
    std::vector<Vector3> positions;
    std::vector<Vector3> pickColors;
    
    unsigned int i = 0;
    unsigned int tris = 0;
    if(mesh->isSimplicial()) {
        tris = mesh->nFaces();
        positions.resize(3*tris);
        pickColors.resize(3*tris);
        
        // Add faces
        for (FacePtr f: mesh->faces()) {
            unsigned int fIdx = fIndices[f];
            HalfedgePtr he = f.halfedge();
            Vector3 color = elementColor(fIdx);
            
            for (unsigned int j = 0; j < 3; j++) {
                VertexPtr v = he.vertex();
                positions[3*i + j] = geometry->position(v);
                pickColors[3*i + j] = color;
                he = he.next();
            }
            i++;
        }
        
    } else {
        for (FacePtr f: mesh->faces()) {
            unsigned int fIdx = fIndices[f];
            HalfedgePtr he = f.halfedge();
            Vector3 color = elementColor(fIdx);
            
            if (f.degree() == 3) {
                for (unsigned int j = 0; j < 3; j++) {
                    VertexPtr v = he.vertex();
                    positions.push_back(geometry->position(v));
                    pickColors.push_back(color);
                    he = he.next();
                }
                i++;
                tris++;
                
            } else {
                std::vector<Triangle> triangles = f.triangulation();
                unsigned int k = (unsigned int)triangles.size();
                for (unsigned int i = 0; i < k; i++ ) {
                    for (unsigned int j = 0; j < 3; j++) {
                        VertexPtr v = triangles[i][j];
                        positions.push_back(geometry->position(v));
                        pickColors.push_back(color);
                    }
                }
                tris += k;
            }
        }
    }
    
    unsigned int verts = 0;
    for (VertexPtr v: mesh->vertices()) {
        verts += v.degree();
        if (v.isBoundary()) verts -= 1;
    }
    
    unsigned int edges = 0;
    for (EdgePtr e: mesh->edges()) {
        edges += 4;
        if (e.isBoundary()) edges -= 2;
    }
    
    size_t size = 3*(tris + edges + verts);
    positions.resize(size);
    pickColors.resize(size);
    
    // Add edges
    for (EdgePtr e: mesh->edges()) {
        unsigned int eIdx = tris + eIndices[e];
        HalfedgePtr h[] = {e.halfedge(), e.halfedge().twin()};
        
        Vector3 n;
        if (h[0].isReal()) n += geometry->normal(h[0].face());
        if (h[1].isReal()) n += geometry->normal(h[1].face());
        n *= (h[0].isReal() && h[1].isReal() ? 0.0005 : 0.001);
        
        Vector3 p1 = geometry->position(h[0].vertex()) + n;
        Vector3 p2 = geometry->position(h[1].vertex()) + n;
        Vector3 color = elementColor(eIdx);
        
        for (int j = 0; j < 2; j++) {
            if (h[j].isReal()) {
                Vector3 p3 = geometry->position(h[j].prev().vertex()) + n;
                
                positions[3*i + 0] = p1;
                positions[3*i + 1] = p2;
                positions[3*i + 2] = p2 + 0.125*(p3 - p2);
                pickColors[3*i + 0] = pickColors[3*i + 1] = pickColors[3*i + 2] = color;
                i++;
                
                positions[3*i + 0] = p1;
                positions[3*i + 1] = positions[3*(i-1) + 2];
                positions[3*i + 2] = p1 + 0.125*(p3 - p1);
                pickColors[3*i + 0] = pickColors[3*i + 1] = pickColors[3*i + 2] = color;
                i++;
            }
        }
    }
    
    // Add vertices
    for (VertexPtr v: mesh->vertices()) {
        unsigned int vIdx = tris + mesh->nEdges() + vIndices[v];
        Vector3 n = 0.0015*geometry->normal(v);
        Vector3 p1 = geometry->position(v) + n;
        Vector3 color = elementColor(vIdx);
        
        for (HalfedgePtr h: v.outgoingHalfedges()) {
            if (h.isReal()) {
                Vector3 p2 = geometry->position(h.twin().vertex()) + n;
                Vector3 p3 = geometry->position(h.prev().vertex()) + n;
                
                positions[3*i + 0] = p1;
                positions[3*i + 1] = p1 + 0.125*(p2 - p1);
                positions[3*i + 2] = p1 + 0.125*(p3 - p1);
                pickColors[3*i + 0] = pickColors[3*i + 1] = pickColors[3*i + 2] = color;
                i++;
            }
        }
    }

    pickGLProgram->setAttribute("a_position", positions, update);
    pickGLProgram->setAttribute("a_color", pickColors, update);
}

void SceneMesh::setColorDefault() {

    // Ensure we have a shader that supports colors
    delete faceGLProgram;
    faceGLProgram = new GLProgram(&SHINY_COLOR_VERT_SHADER, &SHINY_COLOR_FRAG_SHADER, DrawMode::Triangles);
    fillMeshBuffers(false);
    
    std::vector<Vector3> colors(3*mesh->nFaces());
    
    // Walk the faces building an array of colors
    unsigned int i = 0;
    for(FacePtr f : mesh->faces()) {
        HalfedgePtr he = f.halfedge();
        for(unsigned int j = 0; j < 3; j++) {
            VertexPtr v = he.vertex();
            colors[3*i+j] = RGB_TEAL;
            he = he.next();
        }
        i++;
    }
    
    // Store the values in GL buffers
    faceGLProgram->setAttribute("a_color", colors);
}

void SceneMesh::setColorFromScalar(VertexData<double> const &colorData, Colormap const &cm) {

    // Compute automatic bounds based on the ranges of the data
    double minVal = std::numeric_limits<double>::infinity();
    double maxVal = -std::numeric_limits<double>::infinity(); 
    for(VertexPtr v : mesh->vertices()) {
        minVal = std::min<double>(minVal, colorData[v]);
        maxVal = std::max<double>(maxVal, colorData[v]);
    }
   
    // Slightly expand the ranges to avoid nonsense on constant input
    if(minVal == maxVal) {
        double offset = std::max(std::abs(minVal), std::abs(maxVal)) / 1000000.0;  
        offset = std::max(0.00000001, offset);
        minVal -= offset;
        maxVal += offset;
    }

    
    setColorFromScalar(colorData, minVal, maxVal, cm);
}

void SceneMesh::setColorFromScalar(VertexData<double> const &colorData, double minBound, double maxBound, Colormap const &cm) {

    // Ensure we have a shader that supports colors
    delete faceGLProgram;
    faceGLProgram = new GLProgram(&SHINY_COLOR_VERT_SHADER, &SHINY_COLOR_FRAG_SHADER, DrawMode::Triangles);
    fillMeshBuffers(false);


    std::vector<Vector3> colors(3*mesh->nFaces());
 
    double range = maxBound - minBound;

    // Walk the faces building an array of colors
    unsigned int i = 0;
    for(FacePtr f : mesh->faces()) {
        HalfedgePtr he = f.halfedge();
        for(unsigned int j = 0; j < 3; j++) {
            VertexPtr v = he.vertex();

            double val = colorData[v];
            double adjVal = (val - minBound) / range;
            colors[3*i+j] = cm.getValue(adjVal);  

            he = he.next();
        }
        i++;
    }

    // Store the values in GL buffers
    faceGLProgram->setAttribute("a_color", colors);
}

void SceneMesh::setColorFromScalar(FaceData<double> const &colorData, Colormap const &cm) {

    // Compute automatic bounds based on the ranges of the data
    double minVal = std::numeric_limits<double>::infinity();
    double maxVal = -std::numeric_limits<double>::infinity(); 
    for(FacePtr f : mesh->faces()) {
        minVal = std::min<double>(minVal, colorData[f]);
        maxVal = std::max<double>(maxVal, colorData[f]);
    }
   
    // Slightly expand the ranges to avoid nonsense on constant input
    if(minVal == maxVal) {
        double offset = std::max(std::abs(minVal), std::abs(maxVal)) / 1000000.0;  
        offset = std::max(0.00000001, offset);
        minVal -= offset;
        maxVal += offset;
    }

    
    setColorFromScalar(colorData, minVal, maxVal, cm);
}

void SceneMesh::setColorFromScalar(FaceData<double> const &colorData, double minBound, double maxBound, Colormap const &cm) {

    // Ensure we have a shader that supports colors
    delete faceGLProgram;
    faceGLProgram = new GLProgram(&SHINY_COLOR_VERT_SHADER, &SHINY_COLOR_FRAG_SHADER, DrawMode::Triangles);
    fillMeshBuffers(false);


    std::vector<Vector3> colors(3*mesh->nFaces());
 
    double range = maxBound - minBound;

    // Walk the faces building an array of colors
    unsigned int i = 0;
    for(FacePtr f : mesh->faces()) {
        for(unsigned int j = 0; j < 3; j++) {
            double val = colorData[f];
            double adjVal = (val - minBound) / range;
            colors[3*i+j] = cm.getValue(adjVal);  
        }
        i++;
    }

    // Store the values in GL buffers
    faceGLProgram->setAttribute("a_color", colors);
}

void SceneMesh::setColorFromRGB(VertexData<Vector3> const &colorData) {

    // Ensure we have a shader that supports colors
    delete faceGLProgram;
    faceGLProgram = new GLProgram(&SHINY_COLOR_VERT_SHADER, &SHINY_COLOR_FRAG_SHADER, DrawMode::Triangles);
    fillMeshBuffers(false);

    std::vector<Vector3> colors(3*mesh->nFaces());
 
    // Walk the faces building an array of colors
    unsigned int i = 0;
    for(FacePtr f : mesh->faces()) {
        HalfedgePtr he = f.halfedge();
        for(unsigned int j = 0; j < 3; j++) {
            VertexPtr v = he.vertex();
            colors[3*i+j] = clamp(colorData[v], Vector3{0,0,0}, Vector3{1,1,1});
            he = he.next();
        }
        i++;
    }

    // Store the values in GL buffers
    faceGLProgram->setAttribute("a_color", colors);
}

void SceneMesh::setColorFromRGB(FaceData<Vector3> const &colorData)
{
    // Ensure we have a shader that supports colors
    delete faceGLProgram;
    faceGLProgram = new GLProgram(&SHINY_COLOR_VERT_SHADER, &SHINY_COLOR_FRAG_SHADER, DrawMode::Triangles);
    fillMeshBuffers(false);
    
    std::vector<Vector3> colors(3*mesh->nFaces());
    // Walk the faces building an array of colors
    unsigned int i = 0;
    for (FacePtr f : mesh->faces()) {
        HalfedgePtr he = f.halfedge();
        for (unsigned int j = 0; j < 3; j++) {
            colors[3*i+j] = clamp(colorData[f], Vector3{0,0,0}, Vector3{1,1,1});
            he = he.next();
        }
        i++;
    }
    
    // Store the values in GL buffers
    faceGLProgram->setAttribute("a_color", colors);
}

void SceneMesh::setColorFromRGB(CornerData<Vector3> const &colorData)
{
    // Ensure we have a shader that supports colors
    delete faceGLProgram;
    faceGLProgram = new GLProgram(&SHINY_COLOR_VERT_SHADER, &SHINY_COLOR_FRAG_SHADER, DrawMode::Triangles);
    fillVertexBasedMeshBuffers(false);
    
    std::vector<Vector3> colors(18*mesh->nFaces());
    // Walk the faces building an array of colors
    unsigned int i = 0;
    for (FacePtr f: mesh->faces()) {
        unsigned j = 0;
        for (CornerPtr c: f.adjacentCorners()) {
            for (unsigned int k = 0; k < 3; k++) colors[18*i+6*j+k] = clamp(colorData[c],
                                                                            Vector3{0,0,0},
                                                                            Vector3{1,1,1});
            for (unsigned int k = 0; k < 3; k++) colors[18*i+6*j+3+k] = clamp(colorData[c.next()],
                                                                              Vector3{0,0,0},
                                                                              Vector3{1,1,1});
            j++;
        }
        i++;
    }
    
    // Store the values in GL buffers
    faceGLProgram->setAttribute("a_color", colors);
}

void SceneMesh::setColorFromRGB(EdgeData<Vector3> const &colorData)
{
    // Ensure we have a shader that supports colors
    delete faceGLProgram;
    faceGLProgram = new GLProgram(&SHINY_COLOR_VERT_SHADER, &SHINY_COLOR_FRAG_SHADER, DrawMode::Triangles);
    fillEdgeBasedMeshBuffers(false);
    
    std::vector<Vector3> colors(9*mesh->nFaces());
    
    // Walk the faces building an array of colors
    unsigned int i = 0;
    for (FacePtr f: mesh->faces()) {
        unsigned int j = 0;
        for (HalfedgePtr h: f.adjacentHalfedges()) {
            for (unsigned int k = 0; k < 3; k++) colors[9*i+3*j+k] = clamp(colorData[h.edge()],
                                                                           Vector3{0,0,0},
                                                                           Vector3{1,1,1});
            j++;
        }
        i++;
    }
    
    // Store the values in GL buffers
    faceGLProgram->setAttribute("a_color", colors);
}

void SceneMesh::setCheckerBoardTexture(const CornerData<Vector2>& uvs)
{
    delete faceGLProgram;
    faceGLProgram = new GLProgram(&CHECKER_VERT_SHADER, &CHECKER_FRAG_SHADER, DrawMode::Triangles);
    fillMeshBuffers(false);
    
    std::vector<Vector2> texcoords(3*mesh->nFaces());
    std::vector<Vector3> colors(3*mesh->nFaces());
    
    // Walk the faces building an array of texture coordinates and colors
    unsigned int i = 0;
    for (FacePtr f : mesh->faces()) {
        unsigned int j = 0;
        for (CornerPtr c: f.adjacentCorners()) {
            texcoords[3*i+j] = uvs[c]*0.01;
            colors[3*i+j] = RGB_TEAL;
            
            j++;
        }
        i++;
    }
    
    // Store the values in GL buffers
    faceGLProgram->setAttribute("a_color", colors);
    faceGLProgram->setAttribute("a_texcoord", texcoords);
    
    // Set uniform
    float intensity = 1.0;
    faceGLProgram->setUniform("u_intensity", intensity);
}

void SceneMesh::setSelectedTriangle(unsigned int fIdx, const Vector3 color)
{
    std::vector<Vector3> colors(3);
    
    for (unsigned int i = 0; i < 3; i++) {
        colors[i] = color;
    }
    
    // Store the values in GL buffers
    faceGLProgram->setAttribute("a_color", colors, true, 3*fIdx, 3);
}

void SceneMesh :: update( void )
{
    if (bufferMode == VERTEX_BASED_BUFFERS) fillVertexBasedMeshBuffers(true);
    else if (bufferMode == EDGE_BASED_BUFFERS) fillEdgeBasedMeshBuffers(true);
    else fillMeshBuffers(true);
    //fillMeshBuffersDual(true);
    fillMeshBuffersShadow(true);
    fillMeshBuffersPick(true);
}


Geometry<Euclidean>* SceneMesh::getGeometry() {
    return geometry;
}