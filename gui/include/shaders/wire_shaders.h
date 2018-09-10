#pragma once


static const VertShader WIRE_VERT_SHADER =  {
    
    // uniforms
    {
       {"u_viewMatrix", GLData::Matrix44Float},
       {"u_projMatrix", GLData::Matrix44Float},
    },

    // attributes
    {
        {"a_position", GLData::Vector3Float},
        {"a_baryCoord", GLData::Vector3Float},
    },

    // source
    GLSL(150,
      uniform mat4 u_viewMatrix;
      uniform mat4 u_projMatrix;
      in vec3 a_position;
      in vec3 a_baryCoord;
      out vec3 Position;
      out vec3 BaryCoord;

      void main()
      {
          Position = a_position;
          BaryCoord = a_baryCoord;
          gl_Position = u_projMatrix * u_viewMatrix * vec4(Position,1.);
      }
    )
};

// This shader will pass barycentric coordinates to the fragment shader;
// Currently, however, we need to pass in alternative barycentric coordinates
// to hide edges inside the tessellation of non-triangular polygons.
// static const GeomShader WIRE_GEOM_SHADER = {
//     
//     // uniforms
//     {
//         {"u_viewMatrix", GLData::Matrix44Float},
//         {"u_projMatrix", GLData::Matrix44Float},
//     }, 
// 
//     // attributes
//     {
//     },
// 
//     // source
//     GLSL(150,
//         layout(triangles) in;
//         layout(triangle_strip, max_vertices=3) out;
//         uniform mat4 u_viewMatrix;
//         uniform mat4 u_projMatrix;
//         in vec3 vPosition[];
//         out vec3 Position;
//         out vec3 BaryCoord;
//         void main()   {
//             mat4 PV = u_projMatrix * u_viewMatrix;
// 
//             Position = vPosition[0];
//             BaryCoord = vec3(1,0,0);
//             gl_Position = PV * vec4(Position,1.);
//             EmitVertex();
// 
//             Position = vPosition[1];
//             BaryCoord = vec3(0,1,0);
//             gl_Position = PV * vec4(Position,1.);
//             EmitVertex();
// 
//             Position = vPosition[2];
//             BaryCoord = vec3(0,0,1);
//             gl_Position = PV * vec4(Position,1.);
//             EmitVertex();
// 
//             EndPrimitive();
//         }
//     )
// };

static const FragShader WIRE_FRAG_SHADER = {
    
    // uniforms
    {
        {"u_wireframeWeight", GLData::Float},
        {"u_wireColor", GLData::Vector3Float},
    }, 

    // attributes
    {
    },
    
    // textures 
    {
    },
    
    // output location
    "outputF",
    
    // source 
    GLSL(150,
      uniform float u_wireframeWeight;
      uniform vec3 u_wireColor;
      in vec3 Position;
      in vec3 BaryCoord;
      out vec4 outputF;

      // Forward declarations of methods from <shaders/common.h>
      float getEdgeFactor(vec3 UVW);

      void main()
      {
         float edgeFactor = getEdgeFactor(BaryCoord);
         outputF.rgb = u_wireColor;
         outputF.a = u_wireframeWeight * edgeFactor;
      }

    )
};

