#pragma once


// NOTE: You probably dont' want to include this directly... see shaders.h

static const VertShader SIMPLE_VERT_SHADER = {
    // uniforms
    {
        {"u_viewMatrix", GLData::Matrix44Float },
        {"u_projMatrix", GLData::Matrix44Float },
    }, 

    // attributes
    {
        {"a_position", GLData::Vector3Float},
        {"a_color", GLData::Vector3Float}
    },

    // source
    GLSL(150,
      uniform mat4 u_viewMatrix;
      uniform mat4 u_projMatrix;
      in vec3 a_position;
      in vec3 a_color;
      out vec3 Color;

      void main()
      {
          Color = a_color;
          gl_Position = u_projMatrix * u_viewMatrix * vec4(a_position, 1.0) ;
      }
    )
    
};

static const FragShader SIMPLE_FRAG_SHADER = {
    
    // uniforms
    {
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
        in vec3 Color;
        out vec4 outputF;
        void main()
        {
            outputF = vec4(Color, 1.0);
        }
    )
};
