#pragma once

// Interpolates triangles using "phong tessellation"

static const VertShader PHONG_VERT_SHADER =  {
   
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
    },

    // source
    GLSL(410,
      uniform mat4 u_projMatrix;
      uniform mat4 u_viewMatrix;
      in vec3 a_position;
      in vec3 a_normal;
      in vec3 a_color;
      out vec3 Position;
      out vec3 Normal;
      out vec3 Color;

      void main()
      {
          Position = a_position;
          Normal = a_normal;
          Color = a_color;
          gl_Position = u_projMatrix * u_viewMatrix * vec4( a_position, 1. );
      }
    )
};

static const TessShader PHONG_TESS_SHADER =  {
    
    // uniforms
    {},

    // attributes
    {},

    // source
    GLSL(410,
      layout(vertices = 3) out;
      in vec3 Position[]; out vec3 tPosition[];
      in vec3 Normal[];   out vec3 tNormal[];
      in vec3 Color[];   out vec3 tColor[];

      void main()
      {
         tPosition[gl_InvocationID] = Position[gl_InvocationID];
         tNormal[gl_InvocationID] = Normal[gl_InvocationID];
         tColor[gl_InvocationID] = Color[gl_InvocationID];
         if( gl_InvocationID == 0 ) {
            gl_TessLevelInner[0] = 6;
            gl_TessLevelOuter[0] = 6;
            gl_TessLevelOuter[1] = 6;
            gl_TessLevelOuter[2] = 6;
            // TODO make adaptive
         }
      }
    )
};

static const EvalShader PHONG_EVAL_SHADER =  {
    
    // uniforms
    {
        {"u_viewMatrix", GLData::Matrix44Float},
        {"u_projMatrix", GLData::Matrix44Float},
    }, 

    // attributes
    {},

    // source
    GLSL(410,
      layout(triangles, equal_spacing, ccw) in;
      uniform mat4 u_viewMatrix;
      uniform mat4 u_projMatrix;
      in vec3 tPosition[]; out vec3 Position;
      in vec3 tNormal[];   out vec3 Normal;
      in vec3 tColor[];   out vec3 Color;

      void main()
      {
         vec3 p[3];
         p[0] = tPosition[0];
         p[1] = tPosition[1];
         p[2] = tPosition[2];

         vec3 N[3];
         N[0] = tNormal[0];
         N[1] = tNormal[1];
         N[2] = tNormal[2];

         float d[3];
         d[0] = dot( p[0], N[0] );
         d[1] = dot( p[1], N[1] );
         d[2] = dot( p[2], N[2] );

         vec3 q = gl_TessCoord.x * p[0] +
                  gl_TessCoord.y * p[1] +
                  gl_TessCoord.z * p[2] ;

         vec3 r[3];
         r[0] = q + (dot( p[0]-q,N[0]))*N[0];
         r[1] = q + (dot( p[1]-q,N[1]))*N[1];
         r[2] = q + (dot( p[2]-q,N[2]))*N[2];

         Position = gl_TessCoord.x * r[0] +
                    gl_TessCoord.y * r[1] +
                    gl_TessCoord.z * r[2] ;

         Normal = gl_TessCoord.x * tNormal[0] +
                  gl_TessCoord.y * tNormal[1] +
                  gl_TessCoord.z * tNormal[2] ;

         Color = gl_TessCoord.x * tColor[0] +
                 gl_TessCoord.y * tColor[1] +
                 gl_TessCoord.z * tColor[2] ;

         gl_Position = u_projMatrix * u_viewMatrix * vec4(Position, 1.);
      }
    )
};

static const FragShader PHONG_FRAG_SHADER = {
    
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
    GLSL(410,
      uniform vec3 u_eye;
      uniform vec3 u_light;
      in vec3 Position;
      in vec3 Normal;
      in vec3 Color;
      out vec4 outputF;

      float diffuse( vec3 N, vec3 L );
      float specular( vec3 N, vec3 L, vec3 E, float shininess );
      float fresnel( vec3 N, vec3 E );

      void main()
      {
         vec3 N = normalize( Normal );
         vec3 L = normalize( u_light - Position );
         vec3 E = normalize( u_eye - Position );
         vec3 H = normalize( L + E );
         float specAngle = max(dot(H,N),0.);
         float spec = pow(specAngle,30.);
         //vec3 R = 2.*dot(L,N)*N - L;
         vec3 one = vec3( 1., 1., 1. );
         //outputF.rgb = 0.2*Color + 0.8*diffuse(N,L)*Color + 0.1*specular(N,L,E)*one + 0.2*fresnel(N,E)*Color;
         outputF.rgb = 0.2*Color + 0.8*diffuse(N,L)*Color + 0.3*spec*one + 0.2*fresnel(N,E)*Color;
         outputF.rgb = pow(outputF.rgb, vec3(1./.8));
         outputF.a = 1.0;
      }

    )
};

