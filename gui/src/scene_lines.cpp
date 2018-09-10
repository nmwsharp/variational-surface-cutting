#include <scene_lines.h>

#include <iostream>
#include <tuple>
#include <math.h>
 
#include <shaders.h>
#include <shaders/line_shaders.h>

SceneLines::SceneLines(Viewer &parent_)
    :   SceneObject(parent_),
        glProgram(&PASSTHRU_LINE_VERT_SHADER, &LINE_GEOM_SHADER, &SHINY_LINE_FRAG_SHADER, DrawMode::Lines)
{
}

void SceneLines::draw() {

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
    
    glProgram.setUniform("u_lineRadius", lineRadius);

    // Draw
    glProgram.draw();

}


void SceneLines::setLines(std::vector<std::array<Vector3, 2>> const &lines) {

    int nLines = lines.size();
    
    // Copy positions to a buffer
    std::vector<Vector3> pos(nLines*2);
    for(int i = 0; i < nLines; i++) {
        pos[2*i + 0] = lines[i][0];
        pos[2*i + 1] = lines[i][1];
    }
        
    // Constant color
    std::vector<Vector3> colorData(nLines*2);
    std::fill(colorData.begin(), colorData.end(), lineColor);

    // Store data in buffers
    glProgram.setAttribute("a_position", pos);
    glProgram.setAttribute("a_color", colorData);
}


