#pragma once

static const VertShader SHADOW_VERT_SHADER = {
    // uniforms
    {
    }, 

    // attributes
    {
        {"a_position", GLData::Vector3Float},
        {"a_normal", GLData::Vector3Float},
    },

    // source
    GLSL(150,
        in vec3 a_position;
        in vec3 a_normal;
        out vec3 Position;
        out vec3 Normal;
        out vec2 UV;
        void main()
        {
            Position = a_position;
            Normal = a_normal;
            UV = vec2(0,0);
        }
    )
};

static const FragShader SHADOW_FRAG_SHADER = {
    
    // uniforms
    {}, 

    // attributes
    {},
    
    // textures 
    {
    },
    
    // output location
    "outputF",
 
    // source
    GLSL(150,
        out vec4 outputF;

        void main()
        {
           outputF.rgb = vec3( 0., 0., 1. );
           outputF.a = .1;
        }
    )
};

// ===============================================================

// Traditional edge-based silhouette extraction
// Requires triangles_adjacency input

static const GeomShader SHADOW_GEOM_SHADER_TRIADJ = {
    
    // uniforms
    {
        {"u_viewMatrix", GLData::Matrix44Float},
        {"u_projMatrix", GLData::Matrix44Float},
        {"u_light", GLData::Vector3Float},
    }, 

    // attributes
    {
    },

    // source
    GLSL(150,
        layout(triangles_adjacency) in;
        layout(triangle_strip, max_vertices=18) out;
        uniform mat4 u_viewMatrix;
        uniform mat4 u_projMatrix;
        uniform vec3 u_light;
        in vec3 Position[];
        in vec3 Normal[];
        in vec2 UV[];

        void extrudeSide( vec3 v0, vec3 v1 )
        {
            mat4 PV = u_projMatrix * u_viewMatrix;

            gl_Position = PV * vec4( v0,               1.0 ); EmitVertex();
            gl_Position = PV * vec4( v0 - u_light.xyz, 0.0 ); EmitVertex();
            gl_Position = PV * vec4( v1,               1.0 ); EmitVertex();
            gl_Position = PV * vec4( v1 - u_light.xyz, 0.0 ); EmitVertex();
        
            EndPrimitive(); 
        }

        void main()   {
           mat4 PV = u_projMatrix * u_viewMatrix;

           vec3 p[6];
           p[0] = Position[0].xyz;
           p[1] = Position[1].xyz;
           p[2] = Position[2].xyz;
           p[3] = Position[3].xyz;
           p[4] = Position[4].xyz;
           p[5] = Position[5].xyz;

           vec3 e1 = p[2] - p[0];
           vec3 e2 = p[4] - p[2];
           vec3 e3 = p[0] - p[4];
           vec3 e4 = p[1] - p[0];
           vec3 e5 = p[3] - p[2];
           vec3 e6 = p[5] - p[4];

           vec3 N = cross( e1, e2 );
           vec3 L = u_light.xyz - p[0];

           // Handle only back facing triangles
           if (dot(N, L) < 0) {

              N = cross( e4, e1 );

              if (dot(N, L) >= 0) {
                 extrudeSide( p[0], p[2] );
              }

              N = cross( e5, e2 );
              L = u_light.xyz - p[2];

              if (dot(N, L) >= 0) {
                 extrudeSide( p[2], p[4] );
              }

              N = cross( e6, e3 );
              L = u_light.xyz - p[4];

              if (dot(N, L) >= 0) {
                 extrudeSide( p[4], p[0] );
              }

              // render the front cap
              gl_Position = PV * vec4(p[0], 1.0); EmitVertex();
              gl_Position = PV * vec4(p[2], 1.0); EmitVertex();
              gl_Position = PV * vec4(p[4], 1.0); EmitVertex();
              EndPrimitive();

              // render the back cap
              gl_Position = PV * vec4(p[0] - u_light.xyz, 0.0); EmitVertex();
              gl_Position = PV * vec4(p[4] - u_light.xyz, 0.0); EmitVertex();
              gl_Position = PV * vec4(p[2] - u_light.xyz, 0.0); EmitVertex();        
              EndPrimitive();
           }
       }
    )
};

// ===============================================================

// Silhouette extraction from interpolated vertex normals
// Needs only "triangles" input

static const GeomShader SHADOW_GEOM_SHADER_TRI = {
    
    // uniforms
    {
        {"u_viewMatrix", GLData::Matrix44Float},
        {"u_projMatrix", GLData::Matrix44Float},
        {"u_light", GLData::Vector3Float},
    }, 

    // attributes
    {
    },

    // source
    GLSL(150,
        layout(triangles) in;
        layout(triangle_strip, max_vertices=8) out;
        uniform mat4 u_viewMatrix;
        uniform mat4 u_projMatrix;
        uniform vec3 u_light;
        in vec3 Position[];
        in vec3 Normal[];
        in vec2 UV[];

        void extrudeQuad( mat4 PV, vec3 a, vec3 b, vec3 c, vec3 d )
        {
           vec4 a0 = PV * vec4( a, 1. ); vec4 a1 = PV * vec4( a - u_light.xyz, 0. );
           vec4 b0 = PV * vec4( b, 1. ); vec4 b1 = PV * vec4( b - u_light.xyz, 0. );
           vec4 c0 = PV * vec4( c, 1. ); vec4 c1 = PV * vec4( c - u_light.xyz, 0. );
           vec4 d0 = PV * vec4( d, 1. ); vec4 d1 = PV * vec4( d - u_light.xyz, 0. );
        
           gl_Position = a0; EmitVertex();
           gl_Position = b0; EmitVertex();
           gl_Position = c0; EmitVertex();
           gl_Position = d0; EmitVertex();
           gl_Position = c1; EmitVertex();
           gl_Position = d1; EmitVertex();
           gl_Position = a1; EmitVertex();
           gl_Position = b1; EmitVertex();
        
           EndPrimitive(); 
        }

        void extrudeTri( mat4 PV, vec3 a, vec3 b, vec3 c )
        {
           vec4 a0 = PV * vec4( a, 1. ); vec4 a1 = PV * vec4( a - u_light.xyz, 0. );
           vec4 b0 = PV * vec4( b, 1. ); vec4 b1 = PV * vec4( b - u_light.xyz, 0. );
           vec4 c0 = PV * vec4( c, 1. ); vec4 c1 = PV * vec4( c - u_light.xyz, 0. );
        
           gl_Position = a0; EmitVertex();
           gl_Position = b0; EmitVertex();
           gl_Position = c0; EmitVertex();
           gl_Position = b1; EmitVertex();
           gl_Position = c1; EmitVertex();
           gl_Position = a1; EmitVertex();
        
           EndPrimitive(); 
        }

        void emitCaps( mat4 PV, vec3 a, vec3 b, vec3 c )
        {
           vec4 a0 = PV * vec4( a, 1. ); vec4 a1 = PV * vec4( a - u_light.xyz, 0. );
           vec4 b0 = PV * vec4( b, 1. ); vec4 b1 = PV * vec4( b - u_light.xyz, 0. );
           vec4 c0 = PV * vec4( c, 1. ); vec4 c1 = PV * vec4( c - u_light.xyz, 0. );

           gl_Position = a0; EmitVertex();
           gl_Position = b0; EmitVertex();
           gl_Position = c0; EmitVertex();
           EndPrimitive();

           gl_Position = c1; EmitVertex();
           gl_Position = b1; EmitVertex();
           gl_Position = a1; EmitVertex();        
           EndPrimitive();
        }

        void emitSegment( mat4 PV, vec3 a, vec3 b )
        {
           vec4 a0 = PV * vec4( a, 1. ); vec4 a1 = PV * vec4( a - u_light.xyz, 1. );
           vec4 b0 = PV * vec4( b, 1. ); vec4 b1 = PV * vec4( b - u_light.xyz, 1. );

           gl_Position = a0; EmitVertex();
           gl_Position = b0; EmitVertex();
           EndPrimitive();
        }

        void main()   {
           mat4 PV = u_projMatrix * u_viewMatrix;

           float c0 = dot( Normal[0], u_light - Position[0] );
           float c1 = dot( Normal[1], u_light - Position[1] );
           float c2 = dot( Normal[2], u_light - Position[2] );

           vec3 p0 = Position[0] - .005*Normal[0]; // XXX hack
           vec3 p1 = Position[1] - .005*Normal[1];
           vec3 p2 = Position[2] - .005*Normal[2];

           if( c0*c1 < 0. && c0*c2 < 0. )
           {
              // Compute edge-silhouette intersection points
              float s01 = c0/(c0-c1);
              float s02 = c0/(c0-c2);
              vec3 p01 = (1.-s01)*p0 + s01*p1;
              vec3 p02 = (1.-s02)*p0 + s02*p2;

              // Emit shadow volume geometry
              if( c0 > 0. ) extrudeQuad( PV, p1, p2, p01, p02 );
              else extrudeTri( PV, p0, p01, p02 );
           }
           else if( c1*c2 < 0. && c1*c0 < 0. )
           {
              // Compute edge-silhouette intersection points
              float s12 = c1/(c1-c2);
              float s10 = c1/(c1-c0);
              vec3 p12 = (1.-s12)*p1 + s12*p2;
              vec3 p10 = (1.-s10)*p1 + s10*p0;

              // Emit shadow volume geometry
              if( c1 > 0. ) extrudeQuad( PV, p2, p0, p12, p10 );
              else extrudeTri( PV, p1, p12, p10 );
           }
           else if( c2*c0 < 0. && c2*c1 < 0. )
           {
              // Compute edge-silhouette intersection points
              float s20 = c2/(c2-c0);
              float s21 = c2/(c2-c1);
              vec3 p20 = (1.-s20)*p2 + s20*p0;
              vec3 p21 = (1.-s21)*p2 + s21*p1;

              // Emit shadow volume geometry
              if( c2 > 0. ) extrudeQuad( PV, p0, p1, p20, p21 );
              else extrudeTri( PV, p2, p20, p21 );
           }
           else if( c0 < 0. )
           {
              emitCaps( PV, p0, p1, p2 );
           }
      }
   )
};

// ===============================================================

// Silhouette extraction from interpolated vertex normals
// Needs only "triangles" input

static const GeomShader SILHOUETTE_GEOM_SHADER_TRI = {
    
    // uniforms
    {
        {"u_viewMatrix", GLData::Matrix44Float},
        {"u_projMatrix", GLData::Matrix44Float},
        {"u_light", GLData::Vector3Float},
    }, 

    // attributes
    {
    },

    // source
    GLSL(150,
        layout(triangles) in;
        layout(line_strip, max_vertices=2) out;
        uniform mat4 u_viewMatrix;
        uniform mat4 u_projMatrix;
        uniform vec3 u_light;
        in vec3 Position[];
        in vec3 Normal[];
        in vec2 UV[];

        void emitSegment( mat4 PV, vec3 a, vec3 b )
        {
           vec4 a0 = PV * vec4( a, 1. );
           vec4 b0 = PV * vec4( b, 1. );

           gl_Position = a0; EmitVertex();
           gl_Position = b0; EmitVertex();
           EndPrimitive();
        }

        void main()   {
           mat4 PV = u_projMatrix * u_viewMatrix;

           float c0 = dot( Normal[0], u_light - Position[0] );
           float c1 = dot( Normal[1], u_light - Position[1] );
           float c2 = dot( Normal[2], u_light - Position[2] );

           vec3 p0 = Position[0] - .005*Normal[0]; // XXX hack
           vec3 p1 = Position[1] - .005*Normal[1];
           vec3 p2 = Position[2] - .005*Normal[2];

           if( c0*c1 < 0. && c0*c2 < 0. )
           {
              // Compute edge-silhouette intersection points
              float s01 = c0/(c0-c1);
              float s02 = c0/(c0-c2);
              vec3 p01 = (1.-s01)*p0 + s01*p1;
              vec3 p02 = (1.-s02)*p0 + s02*p2;

              // Emit silhouette segment
              emitSegment( PV, p01, p02 );
           }
           else if( c1*c2 < 0. && c1*c0 < 0. )
           {
              // Compute edge-silhouette intersection points
              float s12 = c1/(c1-c2);
              float s10 = c1/(c1-c0);
              vec3 p12 = (1.-s12)*p1 + s12*p2;
              vec3 p10 = (1.-s10)*p1 + s10*p0;

              // Emit silhouette segment
              emitSegment( PV, p12, p10 );
           }
           else if( c2*c0 < 0. && c2*c1 < 0. )
           {
              // Compute edge-silhouette intersection points
              float s20 = c2/(c2-c0);
              float s21 = c2/(c2-c1);
              vec3 p20 = (1.-s20)*p2 + s20*p0;
              vec3 p21 = (1.-s21)*p2 + s21*p1;

              // Emit silhouette segment
              emitSegment( PV, p20, p21 );
           }
      }
   )
};

