#pragma once

#include <viewer.h>

// An abstract base class for something that can be drawn in the scene
// Any subclass must define a method that draws the object; it can
// optionally define methods that are used to render shadows via
// shadow volumes.
class SceneObject {

    public:

        SceneObject(Viewer &parent_);

        virtual ~SceneObject();

        // Hide copy/move constructors that we probably don't want to use
        SceneObject(const SceneObject& other) = delete;
        SceneObject& operator=(const SceneObject& other) = delete;
        SceneObject(SceneObject&& other) = delete;
        SceneObject& operator=(SceneObject&& other) = delete;

        virtual void draw() = 0;
        virtual void drawLight() {}
        virtual void drawShadow() {}
        virtual void drawPick() {};

        bool receiveShadow = true;

        Viewer &parent;

};
