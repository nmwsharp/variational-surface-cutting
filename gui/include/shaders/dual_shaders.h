#pragma once



static const VertShader PASSTHRU_DUAL_VERT_SHADER = {
    // uniforms
    {
    }, 

    // attributes
    {
        {"a_position", GLData::Vector3Float},
    },

    // source
    GLSL(150,
        in vec3 a_position;
        void main()
        {
            gl_Position = vec4(a_position,1.0);
        }
    )
};

static const GeomShader DUAL_GEOM_SHADER = {
    
    // uniforms
    {
        {"u_viewMatrix", GLData::Matrix44Float},
        {"u_projMatrix", GLData::Matrix44Float}
    },

    // attributes
    {
    },

    // source
    GLSL(150,
        layout(lines_adjacency) in;
        layout(line_strip, max_vertices=11) out;
        uniform mat4 u_viewMatrix;
        uniform mat4 u_projMatrix;

        void main()   {
           mat4 PV = u_projMatrix * u_viewMatrix;

           vec4 m = u_viewMatrix * gl_in[0].gl_Position;
           vec4 a = u_viewMatrix * gl_in[1].gl_Position;
           vec4 b = u_viewMatrix * gl_in[2].gl_Position;
           vec4 c = u_viewMatrix * gl_in[3].gl_Position;

           // Offset to avoid Z-fighting (note that alternatives
           // like glPolygonOffset and modifying depth values will
           // disable early Z-cull)
           m.w += .01;
           a.w += .01;
           b.w += .01;
           c.w += .01;

           gl_Position = u_projMatrix*a; EmitVertex();
           gl_Position = u_projMatrix*m; EmitVertex();
           gl_Position = u_projMatrix*b; EmitVertex();
           EndPrimitive();
           gl_Position = u_projMatrix*m; EmitVertex();
           gl_Position = u_projMatrix*c; EmitVertex();
           EndPrimitive();
        }
       )
    };



static const FragShader DUAL_FRAG_SHADER = {
    
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
        out vec4 outputF;

        void main()
        {
           outputF.rgb = vec3( 0., 0., 0. );
           outputF.a = .8;
        }
    )
};

