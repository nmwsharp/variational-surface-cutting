#pragma once

#include <cstdlib>
#include <vector>
#include <memory>

#include <scene_object.h>
#include <gl_utils.h>
#include <viewer.h>
#include <vector3.h>
#include <shaders.h>
#include <geometry.h>


class SceneVectors: public SceneObject {

    public:

        SceneVectors(Viewer &parent);

        virtual void draw() override;
      
        // Options
        bool enabled = true;
        Vector3 arrowColor = RGB_RED;        

        // Set data
        void setVectors(const std::vector<Vector3> &basePoints, const std::vector<Vector3> &vectors, double radius = -1.0);
        void setVectors(Geometry<Euclidean>* geometry, const VertexData<Vector3> &vertexVectors, bool autoscale=false);
        void setVectors(Geometry<Euclidean>* geometry, const FaceData<Vector3> &faceVectors, bool autoscale=false);
       
        // Symmetric vector fields
        void setSymmetricVectors(Geometry<Euclidean>* geometry, const VertexData<Complex> &vertexVectors, int nSym, bool autoscale=false);
        void setSymmetricVectors(Geometry<Euclidean>* geometry, const VertexData<Vector3> &vertexVectors, int nSym, bool autoscale=false);

        // More unusual vector fields
        void setEdgePerpVectors(Geometry<Euclidean>* geometry, const EdgeData<double> &edgeVectors, bool autoscale=false);
        void setOneForm(Geometry<Euclidean>* geometry, const EdgeData<double> &form, bool autoscale=false, bool zeroSparse = false);
        void setEdgeVectors(Geometry<Euclidean>* geometry, const EdgeData<double> &form, bool autoscale=false);

    protected:

        // Drawing options
        GLProgram vectorProgram;

        // Helpers 
        double computeNiceMeshLengthScale(Geometry<Euclidean>* geometry);
};
