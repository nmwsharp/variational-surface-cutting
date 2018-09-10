#pragma once



static const VertShader MULTI_REGION_VERT_SHADER =  {
    
    // uniforms
    {
        {"u_viewMatrix", GLData::Matrix44Float},
        {"u_projMatrix", GLData::Matrix44Float},
    }, 

    // attributes
    {
        {"a_position", GLData::Vector3Float},
        {"a_normal", GLData::Vector3Float},
        {"a_regionVals", GLData::Vector3Float},
        {"a_regionColor0", GLData::Vector3Float},
        {"a_regionColor1", GLData::Vector3Float},
        {"a_regionColor2", GLData::Vector3Float},
    },

    // source
    GLSL(150,
      uniform mat4 u_viewMatrix;
      uniform mat4 u_projMatrix;

      in vec3 a_position;
      in vec3 a_normal;
      in vec3 a_regionVals;
      in vec3 a_regionColor0;
      in vec3 a_regionColor1;
      in vec3 a_regionColor2;

      out vec3 Normal;
      out vec3 Position;
      out vec3 regionVals;
      out vec3 regionColor0;
      out vec3 regionColor1;
      out vec3 regionColor2;

      void main()
      {
          //StripeCoord = a_stripeCoord;
          Position = a_position;
          Normal = a_normal;
          gl_Position = u_projMatrix * u_viewMatrix * vec4(a_position,1.0) ;
          regionVals = a_regionVals;
          regionColor0 = a_regionColor0;
          regionColor1 = a_regionColor1;
          regionColor2 = a_regionColor2;
      }
    )
};

static const FragShader MULTI_REGION_FRAG_SHADER = {
    
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
      uniform int nRegions;

      in vec3 Normal;
      in vec3 Position;
      in vec3 BaryCoord;
      in vec3 regionVals;
      in vec3 regionColor0;
      in vec3 regionColor1;
      in vec3 regionColor2;

      out vec4 outputF;

      // Forward declarations of methods from <shaders/common.h>
      vec3 getSurfaceColor();
      vec4 lightSurface( vec3 position, vec3 normal, vec3 color, vec3 light, vec3 eye );

      vec3 minRegionColor() {

        if(regionVals[0] < regionVals[1] && regionVals[0] < regionVals[2]) {
            return regionColor0;
        }
        if(regionVals[1] < regionVals[2]) {
            return regionColor1;
        }
        return regionColor2;
      }

      void main()
      {

         vec3 Color = minRegionColor();
         outputF = lightSurface( Position, Normal, Color, u_light, u_eye );
      }

    )
};

static const VertShader TWO_REGION_VERT_SHADER =  {
    
    // uniforms
    {
        {"u_viewMatrix", GLData::Matrix44Float},
        {"u_projMatrix", GLData::Matrix44Float},
    }, 

    // attributes
    {
        {"a_position", GLData::Vector3Float},
        {"a_normal", GLData::Vector3Float},
        {"a_region", GLData::Float},
        {"a_colorPos", GLData::Vector3Float},
        {"a_colorNeg", GLData::Vector3Float},
    },

    // source
    GLSL(150,
      uniform mat4 u_viewMatrix;
      uniform mat4 u_projMatrix;

      in vec3 a_position;
      in vec3 a_normal;
      in float a_region;
      in vec3 a_colorPos;
      in vec3 a_colorNeg;

      out vec3 Normal;
      out vec3 Position;
      out float region;
      out vec3 colorPos;
      out vec3 colorNeg;

      void main()
      {
          //StripeCoord = a_stripeCoord;
          Position = a_position;
          Normal = a_normal;
          gl_Position = u_projMatrix * u_viewMatrix * vec4(a_position,1.0) ;
          region = a_region;
          colorPos = a_colorPos;
          colorNeg = a_colorNeg;
          
      }
    )
};

static const FragShader TWO_REGION_FRAG_SHADER = {
    
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

      in vec3 Normal;
      in vec3 Position;
      in vec3 BaryCoord;
      in float region;
      in vec3 colorPos;
      in vec3 colorNeg;

      out vec4 outputF;

      // Forward declarations of methods from <shaders/common.h>
      vec3 getSurfaceColor();
      vec4 lightSurface( vec3 position, vec3 normal, vec3 color, vec3 light, vec3 eye );

      void main()
      {

         vec3 Color;
         if(region > 0) {
             Color = colorPos;
         } else {
             Color = colorNeg;
         }

         outputF = lightSurface( Position, Normal, Color, u_light, u_eye );
      }

    )
};



// NOTE: You probably dont' want to include this directly... see shaders.h
static const VertShader CUT_CHECKER_VERT_SHADER =  {
    
    // uniforms
    {
        {"u_viewMatrix", GLData::Matrix44Float},
        {"u_projMatrix", GLData::Matrix44Float},
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

static const FragShader CUT_CHECKER_FRAG_SHADER = {
    
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
      vec4 lightSurface( vec3 position, vec3 normal, vec3 color, vec3 light, vec3 eye );
      float check( vec2 uv );

      void main()
      {
         float intensity = check(Texcoord) + .25;
         vec3 checkColor = Color * intensity;
         outputF = lightSurface( Position, Normal, checkColor, u_light, u_eye );
      }

    )
};

