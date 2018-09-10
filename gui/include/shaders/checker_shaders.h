#pragma once


// NOTE: You probably dont' want to include this directly... see shaders.h
static const VertShader CHECKER_VERT_SHADER =  {
    
    // uniforms
    {
        {"u_viewMatrix", GLData::Matrix44Float},
        {"u_projMatrix", GLData::Matrix44Float},
        {"u_intensity", GLData::Float},
    }, 

    // attributes
    {
        {"a_position", GLData::Vector3Float},
        {"a_color", GLData::Vector3Float},
        {"a_normal", GLData::Vector3Float},
        {"a_texcoord", GLData::Vector2Float},
    },

    // source
    GLSL(150,
      uniform mat4 u_viewMatrix;
      uniform mat4 u_projMatrix;
      in vec3 a_position;
      in vec3 a_normal;
      in vec4 a_color;
      in vec2 a_texcoord;
      out vec3 Color;
      out vec3 Normal;
      out vec3 Position;
      out vec2 Texcoord;

      void main()
      {
          Color = a_color.xyz;
          Position = a_position;
          Normal = a_normal;
          Texcoord = a_texcoord;
          gl_Position = u_projMatrix * u_viewMatrix * vec4(a_position,1.0) ;
      }
    )
};

static const FragShader CHECKER_FRAG_SHADER = {
    
    // uniforms
    {
        {"u_eye", GLData::Vector3Float},
        {"u_light", GLData::Vector3Float},
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
      uniform vec3 u_eye;
      uniform vec3 u_light;
      uniform float u_intensity;
      in vec3 Position;
      in vec3 Normal;
      in vec3 Color;
      in vec2 Texcoord;
      out vec4 outputF;

      // Declare methods from <shaders/common.h>
      float check( vec2 uv );
      vec4 lightSurface( vec3 position, vec3 normal, vec3 color, vec3 light, vec3 eye );

      void main()
      {
         float d = sqrt(25.6*length( Texcoord ));
         vec3 checkColor = clamp(1.5-d,0,1.5) * Color * check(Texcoord) * u_intensity;
         outputF = lightSurface( Position, Normal, checkColor, u_light, u_eye );
      }

    )
};

