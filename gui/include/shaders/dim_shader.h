#pragma once

// This shader draws a full-screen black quad with opacity specified by the float parameter u_dimAmount.
// It ignores any input geometry, using the geometry shader to generate two output triangles.

static const VertShader DIM_VERT_SHADER = {
    // uniforms
    {}, 

    // attributes
    {
       {"a_position", GLData::Vector3Float},
    },

    // source
    GLSL(150,
      void main() {}
    )
    
};

static const GeomShader DIM_GEOM_SHADER = {
    
    // uniforms
    {}, 

    // attributes
    {},

    // source
    GLSL(150,
        layout(points) in;
        layout(triangle_strip, max_vertices=4) out;
        void main() {
           gl_Position = vec4( -1., -1., 1., 1. ); EmitVertex();
           gl_Position = vec4(  1., -1., 1., 1. ); EmitVertex();
           gl_Position = vec4( -1.,  1., 1., 1. ); EmitVertex();
           gl_Position = vec4(  1.,  1., 1., 1. ); EmitVertex();
           EndPrimitive();
        }
    )
};

static const FragShader DIM_FRAG_SHADER = {
    
    // uniforms
    {
        {"u_dimAmount", GLData::Float },
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
        uniform float u_dimAmount;
        void main()
        {
            outputF = vec4(0., 0., 0., u_dimAmount);
        }
    )
};

