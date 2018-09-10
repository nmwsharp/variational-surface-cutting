#pragma once


// NOTE: You probably don't want to include this directly... see shaders.h


static const VertShader PASSTHRU_LINE_VERT_SHADER = {
    // uniforms
    {
    }, 

    // attributes
    {
        {"a_position", GLData::Vector3Float},
        {"a_color", GLData::Vector3Float},
    },

    // source
    GLSL(150,
        in vec3 a_position;
        in vec3 a_color;
        out vec3 Color;
        void main()
        {
            Color = a_color;
            gl_Position = vec4(a_position,1.0);
        }
    )
};


static const GeomShader LINE_GEOM_SHADER = {
    
    // uniforms
    {
        {"u_viewMatrix", GLData::Matrix44Float},
        {"u_projMatrix", GLData::Matrix44Float},
        {"u_lineRadius", GLData::Float},
    }, 

    // attributes
    {
    },

    // source
    GLSL(150,
        layout(lines) in;
        layout(triangle_strip, max_vertices=100) out;
        in vec3 Color[];
        uniform mat4 u_viewMatrix;
        uniform mat4 u_projMatrix;
        uniform float u_lineRadius;
        out vec3 colorToFrag;
        out vec3 worldNormalToFrag;
        out vec3 worldPosToFrag;
        void main()   {
            mat4 PV = u_projMatrix * u_viewMatrix;
            const int nTheta = 8;
            const float PI = 3.14159265358;
            const float delTheta = 2*PI / nTheta;
            const vec3 anyDir = vec3(0.1293832, -0.876434, 0.55236); // provably random

            vec3 pos0 = gl_in[0].gl_Position.xyz;
            vec3 pos1 = gl_in[1].gl_Position.xyz;
            vec3 dir = normalize(pos1 - pos0);

            vec3 yBasis = normalize(cross(dir, anyDir));
            vec3 xBasis = cross(dir, yBasis);

            for (int iTheta = 0; iTheta <= nTheta; iTheta++) { /* duplicate first/last point to complete strip */

                float theta = delTheta * iTheta;
                float cosTheta = cos(theta);
                float sinTheta = sin(theta);

                vec3 normal = cosTheta * xBasis + sinTheta * yBasis;

                // low point
                vec3 low = pos0 + normal*u_lineRadius;
                vec4 worldPosLower = vec4(low, 1);
                gl_Position = PV * worldPosLower;
                worldPosToFrag = worldPosLower.xyz;
                worldNormalToFrag = normal;
                colorToFrag = Color[0];
                EmitVertex();

                // high point
                vec3 high = pos1 + normal*u_lineRadius;
                vec4 worldPosUpper = vec4(high, 1);
                gl_Position = PV * worldPosUpper;
                worldPosToFrag = worldPosUpper.xyz;
                worldNormalToFrag = normal;
                colorToFrag = Color[0];
                EmitVertex();

            }

            EndPrimitive();
        }

    )
};



static const FragShader SHINY_LINE_FRAG_SHADER = {
    
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
        in vec3 colorToFrag;
        in vec3 worldNormalToFrag;
        in vec3 worldPosToFrag;
        out vec4 outputF;

        // Forward declarations of methods from <shaders/common.h>
        vec4 lightSurface( vec3 position, vec3 normal, vec3 color, vec3 light, vec3 eye );

        void main()
        {
           outputF = lightSurface(worldPosToFrag, worldNormalToFrag, colorToFrag, u_light, u_eye );
        }
    )
};

