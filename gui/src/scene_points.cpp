#include <scene_points.h>

#include <iostream>
#include <tuple>
#include <math.h>
 
#include <shaders.h>
#include <shaders/sphere_shaders.h>

ScenePoints::ScenePoints(Viewer &parent_)
    :   SceneObject(parent_),
        glProgram(&PASSTHRU_SPHERE_VERT_SHADER, &SPHERE_GEOM_SHADER, &SHINY_SPHERE_FRAG_SHADER, DrawMode::Points)
{
}

void ScenePoints::draw() {

    if(!enabled) {
        return;
    }

    // Set uniforms
    glm::mat4 viewMat = parent.camera.getViewMatrix();
    glProgram.setUniform("u_viewMatrix", glm::value_ptr(viewMat));

    glm::mat4 projMat = parent.camera.getPerspectiveMatrix();
    glProgram.setUniform("u_projMatrix", glm::value_ptr(projMat));

    Vector3 eyePos = parent.camera.getCameraWorldPosition();
    glProgram.setUniform("u_eye", eyePos);

    Vector3 lightPos = parent.camera.getLightWorldPosition();
    glProgram.setUniform("u_light", lightPos);

    // Draw
    glProgram.draw();

}

// Compute a sphere of radius 0 centered at the origin. Returns (position, normal, index).
void canonicalSphere(int nTheta, int nPhi, float *pos, float *norm, unsigned int* ind) {

    // TODO This should be refactored to use Vector3's rather than buffers of floats

    float delPhi = PI / (nPhi + 1);
    float delTheta = 2*PI / nTheta;

    // All the points in the middle
    for(int iPhi = 0; iPhi < nPhi; iPhi++) {

        float phi = -PI/2.0 + (iPhi+1) * delPhi;

        for(int iTheta = 0; iTheta < nTheta; iTheta++) {

            float theta = delTheta * iTheta;
            int thisPtNum = iPhi * nTheta + iTheta;

            // Compute positions
            float x = cos(phi) * cos(theta);
            float y = cos(phi) * sin(theta);
            float z = sin(phi);
            pos[3*(thisPtNum) + 0] = x;
            pos[3*(thisPtNum) + 1] = y;
            pos[3*(thisPtNum) + 2] = z;

            // Fill the index array to make little triangles out of each square
            if(iPhi != (nPhi-1)) {
                int offSet = 3*(2*nTheta*iPhi + 2*iTheta);

                int nextTheta = (iTheta == (nTheta-1)) ? (-nTheta + 1) : 1;

                // First triangle
                ind[offSet + 0] = thisPtNum;
                ind[offSet + 1] = thisPtNum + nTheta;
                ind[offSet + 2] = thisPtNum + nTheta + nextTheta;

                // Second triangle
                ind[offSet + 3] = thisPtNum;
                ind[offSet + 4] = thisPtNum + nTheta + nextTheta;
                ind[offSet + 5] = thisPtNum + nextTheta;
            }

        }
    }


    // Special points as the north and south pole
    pos[3*(nTheta*nPhi) + 0] = 0.0;
    pos[3*(nTheta*nPhi) + 1] = 0.0;
    pos[3*(nTheta*nPhi) + 2] = -1.0; // south pole

    pos[3*(nTheta*nPhi + 1) + 0] = 0.0;
    pos[3*(nTheta*nPhi + 1) + 1] = 0.0;
    pos[3*(nTheta*nPhi + 1) + 2] = 1.0; // north pole


    // Index triangles to the north and south pole
    int offset = 3*(2*nTheta*(nPhi-1));
    int lastPoint = nTheta * nPhi;
    for(int iTheta = 0; iTheta < nTheta; iTheta++) {

        int nextTheta = (iTheta == (nTheta-1)) ? (-nTheta + 1) : 1;

        // To south pole
        ind[offset + 3*iTheta + 0] = nTheta*nPhi;
        ind[offset + 3*iTheta + 1] = iTheta;
        ind[offset + 3*iTheta + 2] = iTheta + nextTheta;

        // To north pole
        ind[offset + 3*(nTheta + iTheta) + 0] = nTheta*nPhi + 1;
        ind[offset + 3*(nTheta + iTheta) + 1] = nTheta*(nPhi-1) + iTheta;
        ind[offset + 3*(nTheta + iTheta) + 2] = nTheta*(nPhi-1) + iTheta + nextTheta;
    }

    // Normals are conveniently identical to positions
    std::copy(pos, pos + 3*(nTheta*nPhi+2), norm);
}

void ScenePoints::setPositions(std::vector<Vector3> const &pos) {

    // Constant color
    std::vector<Vector3> colorData(pos.size());
    std::fill(colorData.begin(), colorData.end(), pointColor);
   
    // Constant size
    std::vector<double> sizeData(pos.size());
    std::fill(sizeData.begin(), sizeData.end(), sphereRad);

    // Store data in buffers
    glProgram.setAttribute("a_position", pos);
    glProgram.setAttribute("a_color", colorData);
    glProgram.setAttribute("a_pointRadius", sizeData);

}


