#version 410 core
layout (location = 0) in vec3 inPosition;

layout (std140) uniform Transform {
    mat4 projection;
    mat4 view;
    mat4 model;
};

void main()
{
    gl_Position = projection * view * model * vec4(inPosition, 1.0);
}
