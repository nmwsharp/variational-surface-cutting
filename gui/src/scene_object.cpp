#include <scene_object.h>


SceneObject::SceneObject(Viewer &parent_)
    : parent(parent_)
{
    parent.addSceneObject(this);
}

SceneObject::~SceneObject() {
    parent.removeSceneObject(this);
}

