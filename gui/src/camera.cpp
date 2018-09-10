#include <camera.h>

#include <iostream>
#include <cmath>

#include <utilities.h>

Camera::Camera() {
}

void Camera::mouseDragEvent(double oldX, double oldY, double newX, double newY, bool isRotating) {
  
    if(oldX == newX && oldY == newY) {
        return;
    }
 
    // Convert to [-1,1] screen coordinates
    oldX = oldX*2 - 1.0;
    oldY = -(oldY*2 - 1.0);
    newX = newX*2 - 1.0;
    newY = -(newY*2 - 1.0);

    double delX = newX - oldX;
    double delY = newY - oldY;

    // Process a rotation
    if(isRotating) {

        double delTheta = -delX * PI;
        double delPhi = delY * PI;

        Vector3 leftDirection = cross(upDirection, cameraDirection);

        // Rotate corresponding to theta
        // (the 'up' direction would never change, since it's what we're rotating around)
        cameraDirection = cameraDirection.rotate_around(upDirection, delTheta);

        // Rotate corresponding to phi
        cameraDirection = cameraDirection.rotate_around(leftDirection, delPhi);
        upDirection = upDirection.rotate_around(leftDirection, delPhi);

        /* Switch to the below method for a nifty arcball-like effect. It works, but never quite felt right.
        
        // When clicking near the center following a theta-phi rotation policy. When clicking near the outside, follow twisting policy
        double clickRad = clamp(0.9*sqrt(newX*newX + newY*newY), 0.0, 1.0);
        double rotationBlendFactor = clamp(clickRad, 0.0, 1.0);

        // Adjust angles for the theta-phi policy
        double delTheta = -delX * PI * (1.0-rotationBlendFactor);
        double delPhi = delY * PI * (1.0-rotationBlendFactor);

        Vector3 leftDirection = cross(upDirection, cameraDirection);

        // Rotate corresponding to theta
        // (the 'up' direction would never change, since it's what we're rotating around)
        cameraDirection = cameraDirection.rotate_around(upDirection, delTheta);

        // Rotate corresponding to phi
        cameraDirection = cameraDirection.rotate_around(leftDirection, delPhi);
        upDirection = upDirection.rotate_around(leftDirection, delPhi);


    
        // Rotate corresponding to a twisting motion to mimic the nice part of arcball
        double rotateAmount = -asin(cross(unit(Vector3{newX, newY, 0}), unit(Vector3{delX, delY, 0})).z) * rotationBlendFactor / 50;
        // The normalizations above can lead to invalid values. Sidestep the problem the lazy way.
        if(std::isfinite(rotateAmount)) {
            upDirection = upDirection.rotate_around(cameraDirection, rotateAmount);
        }
        */
    } 
    // Process a translation
    else {
        double movementScale = dataLengthScale * 0.3;
        cameraSpaceTranslate -= movementScale*Vector3{delX, delY, 0}; 
    }

}

void Camera::mouseScrollEvent(double scrollAmount, bool scrollClipPlane) {

    // Adjust the near clipping plane
    if(scrollClipPlane) {
        if(scrollAmount > 0) {
            nearClipRatio -= 0.03 * nearClipRatio; 
        } else {
            nearClipRatio += 0.03 * nearClipRatio; 
        }

    } 
    // Translate the camera forwards and backwards 
    else {
        if(scrollAmount > 0) {
            dist -= 0.03 * dataLengthScale; 
        } else {
            dist += 0.03 * dataLengthScale; 
        }
    }
    
}

void Camera::setWindowSize(int width_, int height_) {
    width = width_;
    height = height_;
    
}


glm::mat4 Camera::getViewMatrix() const {


    // Map from world coordinates to camera coordinates
    glm::vec3 eye(cameraDirection.x, cameraDirection.y, cameraDirection.z);
    glm::vec3 center (lookAtPoint.x, lookAtPoint.y, lookAtPoint.z);
    glm::vec3 up(upDirection.x, upDirection.y, upDirection.z);
    glm::mat4 rotateToCameraCoordinates = glm::lookAt(center + eye, center, up);
    
   
    // Translate in camera space
    // This is used to support the "sliding" in the screen plane by dragging while holding 'shift'
	glm::mat4 translateCamera = glm::translate(glm::mat4(1.0), glm::vec3(-cameraSpaceTranslate.x, -cameraSpaceTranslate.y, -dist - cameraSpaceTranslate.z));

    // Compose the operations together
	glm::mat4 view = translateCamera * rotateToCameraCoordinates;

    return view;
}



glm::mat4 Camera::getPerspectiveMatrix() const {

    //double farClip = farClipRatio * dist;
    double farClip = farClipRatio * dataLengthScale;
    double nearClip = nearClipRatio * dataLengthScale;
	double fovRad = glm::radians(fov);
    double aspectRatio = (float) width / height;
	
	return glm::perspective(fovRad, aspectRatio, nearClip, farClip);
}

glm::mat4 Camera::getOrthoMatrix() const {
 
    double farClip = farClipRatio * dataLengthScale;
    double nearClip = nearClipRatio * dataLengthScale;
    
    return glm::ortho(-1.0, 1.0, -1.0, 1.0, nearClip, farClip);
}

Vector3 Camera::getCameraWorldPosition() const {
    //return lookAtPoint + cameraDirection * dist;

   // This will work no matter how the view matrix is constructed...
    glm::mat4 invViewMat = inverse( getViewMatrix() );
    return Vector3{ invViewMat[3][0], invViewMat[3][1], invViewMat[3][2] };
}

Vector3 Camera::getLightWorldPosition() const {
    // The light is over the right shoulder of the camera
    Vector3 leftHand = unit(cross(cameraDirection, upDirection));
    Vector3 lightPos = lookAtPoint + 10 * (dist + dataLengthScale) * cameraDirection;
    lightPos = lightPos.rotate_around(upDirection, 3.14159 / 8);
    lightPos = lightPos.rotate_around(leftHand, 3.14159 / 8);

    // Add a vertical component to the light so that when the
    // camera is at the horizon, the shadow doesn't shoot off to infinity
    lightPos = ( lightPos + Vector3{0.,20.,0.} ) / 2.;

    return lightPos;
}
