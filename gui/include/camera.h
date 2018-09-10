#pragma once

#include <vector3.h>

#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtc/matrix_transform.hpp> // glm::translate, glm::rotate, glm::scale, glm::perspective
#include <glm/gtc/constants.hpp> // glm::pi
#include <glm/gtc/type_ptr.hpp> // glm::value_ptr
#include <glm/gtx/string_cast.hpp> // glm::to_string


class Camera {

    public:
    
        Camera(); 
        glm::mat4 getViewMatrix() const;
        glm::mat4 getPerspectiveMatrix() const;
        glm::mat4 getOrthoMatrix() const;
        Vector3 getCameraWorldPosition() const;
        Vector3 getLightWorldPosition() const;

        void mouseDragEvent(double oldX, double oldY, double newX, double newY, bool isRotating);
        void mouseScrollEvent(double scrollAmount, bool scrollClipPlane);
        void setWindowSize(int width, int height);
        
        double dist = 1.0;
        double dataLengthScale = 1.0; // a representative scale for the data, independent of the current view distance
        Vector3 lookAtPoint {0,0,0};
        Vector3 cameraSpaceTranslate {0,0,0};
        Vector3 cameraDirection {0,0,1};
        Vector3 upDirection {0,1,0};

    private:
        
        int width = -1;
        int height = -1;
        double fov = 65.0;
        double nearClipRatio = 0.1;
        double farClipRatio = 20.0;
};
