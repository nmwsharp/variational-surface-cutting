#pragma once

#include <cstdlib>
#include <vector>
#include <memory>

#include <scene_object.h>
#include <gl_utils.h>
#include <colormaps.h>
#include <viewer.h>
#include <vector3.h>
#include <shaders.h>
#include <geometry.h>

#define VERTEX_BASED_BUFFERS 0
#define EDGE_BASED_BUFFERS 1
#define FACE_BASED_BUFFERS 2

class SceneMesh: public SceneObject {

    public:

        SceneMesh(Viewer &parent, Geometry<Euclidean>* geometryToDraw);
        ~SceneMesh();
        std::string name = "unnamed geometry";

        virtual void draw() override;
        virtual void drawLight() override;
        virtual void drawShadow() override;
        virtual void drawPick() override;
        
        // Coloring-related
        void setColorDefault();
        void setColorFromScalar(VertexData<double> const &colorData, Colormap const &cm = CM_VIRIDIS);
        void setColorFromScalar(VertexData<double> const &colorData, double minBound, double maxBound, Colormap const &cm = CM_VIRIDIS);
        void setColorFromScalar(FaceData<double> const &colorData, Colormap const &cm = CM_VIRIDIS);
        void setColorFromScalar(FaceData<double> const &colorData, double minBound, double maxBound, Colormap const &cm = CM_VIRIDIS);
        void setColorFromRGB(VertexData<Vector3> const &colorData);
        void setColorFromRGB(FaceData<Vector3> const &colorData);
        void setColorFromRGB(CornerData<Vector3> const &colorData);
        void setColorFromRGB(EdgeData<Vector3> const &colorData);
    
        // Texture-related
        void setCheckerBoardTexture(const CornerData<Vector2>& uvs);
    
        // Selection-related
        void setSelectedTriangle(unsigned int fIdx, const Vector3 color);
    
        // Adjust the camera for this data
        void resetCameraForData();

        // Should be called whenever the mesh is modified
        void update( void );

        Geometry<Euclidean>* getGeometry();

    protected:

        // Drawing options
        GLProgram *faceGLProgram, *wireGLProgram, *dualGLProgram,
                  *lightGLProgram, *shadowGLProgram, *pickGLProgram;
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geometry;
        Vector3 dataCenter;
        double dataScale;
        int bufferMode;

        void updateGlBuffers(std::vector<Vector3> positions, std::vector<Vector3> normals,
                             std::vector<Vector3> vertexNormals, std::vector<Vector3> barycentric,
                             bool update);
        virtual void fillMeshBuffers(bool update = false);
        virtual void fillVertexBasedMeshBuffers(bool update = false);
        virtual void fillEdgeBasedMeshBuffers(bool update = false);
        virtual void fillMeshBuffersDual(bool update = false);
        virtual void fillMeshBuffersShadow(bool update = false);
        virtual void fillMeshBuffersPick(bool update = false);
};
