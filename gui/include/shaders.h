#pragma once

#include <cstdlib>
#include <string>
#include <vector>

// Make syntax  nicer like this, but we lose line numbers in GL debug output
#define GLSL(version, shader)  "#version " #version "\n" #shader


// Enum for openGL data types
enum class GLData {Vector2Float, Vector3Float, Vector4Float, Matrix44Float, Float, Int, UInt, Index};

// Types to encapsulate uniform and attribute variables in shaders
struct ShaderUniform {
    const std::string name;
    const GLData type;
};
struct ShaderAttribute {
    const std::string name;
    const GLData type;
};
struct ShaderTexture {
    const std::string name;
    const unsigned int dim;
};


// Types which represents shaders and the values they require
struct VertShader {
    const std::vector<ShaderUniform> uniforms;
    const std::vector<ShaderAttribute> attributes;
    const std::string src;
};
struct TessShader {
    const std::vector<ShaderUniform> uniforms;
    const std::vector<ShaderAttribute> attributes;
    const std::string src;
};
struct EvalShader {
    const std::vector<ShaderUniform> uniforms;
    const std::vector<ShaderAttribute> attributes;
    const std::string src;
};
struct GeomShader {
    const std::vector<ShaderUniform> uniforms;
    const std::vector<ShaderAttribute> attributes;
    const std::string src;
};
struct FragShader {
    const std::vector<ShaderUniform> uniforms;
    const std::vector<ShaderAttribute> attributes;
    const std::vector<ShaderTexture> textures;
    const std::string outputLoc;
    const std::string src;
};


// === Welcome to Nick's Shader Emporium ===

/* 
 * = Simple shaders:
 *    - (vert) SIMPLE_VERT_SHADER
 *    - (frag) SIMPLE_FRAG_SHADER
 * 
 * = Shiny shaders:
 *    - (vert) SHINY_VERT_SHADER 
 *    - (frag) SHINY_FRAG_SHADER
 *
 * = Stripe shaders:
 *    - (vert) STRIPE_SHINY_VERT_SHADER 
 *    - (frag) STRIPE_SHINY_FRAG_SHADER 
 *
 * = Checkerboard shaders:
 *    - (vert) CHECKER_VERT_SHADER 
 *    - (frag) CHECKER_FRAG_SHADER 
 *
 * = Sphere shaders:
 *    - (vert) PASSTHRU_SPHERE_VERT_SHADER
 *    - (geom) SPHERE_GEOM_SHADER
 *    - (frag) SHINY_SPHERE_FRAG_SHADER 
 *
 *
 *
 */

