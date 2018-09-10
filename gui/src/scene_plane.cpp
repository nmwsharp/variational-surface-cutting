#include <scene_plane.h>

#include <iostream>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <limits>
 
#include <shaders.h>
#include <shaders/checker_shaders.h>

ScenePlane::ScenePlane(Viewer &parent_)
    :   SceneObject(parent_),
        showPlane( true ),
        shader(&CHECKER_VERT_SHADER, &CHECKER_FRAG_SHADER, DrawMode::Triangles)
{
   setGeometry( 0., Vector3{ 0., 0., 1. }, Vector3{ 1., 0., 0. }, Vector3 {0,0,0} );
}


void ScenePlane::drawPlane() {

   if( !showPlane )
   {
      return;
   }

    // Set uniforms
    glm::mat4 viewMat = parent.camera.getViewMatrix();
    shader.setUniform("u_viewMatrix", glm::value_ptr(viewMat));

    glm::mat4 projMat = parent.camera.getPerspectiveMatrix();
    shader.setUniform("u_projMatrix", glm::value_ptr(projMat));

    Vector3 eyePos = parent.camera.getCameraWorldPosition();
    shader.setUniform("u_eye", eyePos);

    Vector3 lightPos = parent.camera.getLightWorldPosition();
    shader.setUniform("u_light", lightPos);

    shader.setUniform("u_intensity", intensity);

    // Draw
    shader.draw();
}

void ScenePlane::draw() {
   intensity = 1.;
   drawPlane();
}

void ScenePlane::drawLight() {
   intensity = .5;
   drawPlane();
}

void ScenePlane::setGeometry( double offset, Vector3 e1, Vector3 e2, Vector3 center, double _scale )
{
   scale = _scale;

   // Two triangles make a quad
   int nFaces = 2;
   std::vector<Vector3> positions(3*nFaces);
   std::vector<Vector3> normals(3*nFaces);
   std::vector<Vector3> colors(3*nFaces);
   std::vector<Vector2> texcoords(3*nFaces);

   // Orthonormalize basis
   e1 = unit( e1 );
   e2 -= dot( e2, e1 )*e1;
   e2 = unit( e2 );
   Vector3 e3 = cross( e1, e2 );

   const double S = 1000 * scale;
   positions[0] = center + S * ( -e1 - e2 ) + offset * e3;
   positions[1] = center + S * (  e1 - e2 ) + offset * e3;
   positions[2] = center + S * (  e1 + e2 ) + offset * e3;

   positions[3] = center + S * ( -e1 - e2 ) + offset * e3;
   positions[4] = center + S * (  e1 + e2 ) + offset * e3;
   positions[5] = center + S * ( -e1 + e2 ) + offset * e3;

   texcoords[0] = Vector2{ -1., -1. };
   texcoords[1] = Vector2{  1., -1. };
   texcoords[2] = Vector2{  1.,  1. };

   texcoords[3] = Vector2{ -1., -1. };
   texcoords[4] = Vector2{  1.,  1. };
   texcoords[5] = Vector2{  1., -1. };

   for( unsigned int i = 0; i < 6; i++ )
   {
      // All normals point up
      normals[i] = e3;

      // All vertex colors are white
      colors[i] = Vector3{ 1., 1., 1. };
   }

   // Store the values in GL buffers
   shader.setAttribute("a_position", positions);
   shader.setAttribute("a_normal", normals);
   shader.setAttribute("a_color", colors);
   shader.setAttribute("a_texcoord", texcoords);
}
