#include "RenderData.h"

GLMesh::GLMesh(Mesh& mesh0):
mesh(mesh0),
vao(0),
vbo(0)
{
    
}

void GLMesh::fillBuffers(const std::vector<Eigen::Vector3f>& colors)
{
    vertices.resize(3*(mesh.faces.size() - mesh.boundaries.size()));
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int i = 0;
            int fIdx = 3*f->index;
            
            HalfEdgeCIter h = f->he;
            do {
                VertexCIter v = h->vertex;
                
                // set vertex
                vertices[fIdx+i].position = v->position.cast<float>();
                vertices[fIdx+i].normal = v->normal().cast<float>();
                vertices[fIdx+i].color = colors[v->index];
                vertices[fIdx+i].uv = v->uv.cast<float>();
                i++;
                
                h = h->next;
            } while (h != f->he);
            
            // set barycenters
            vertices[fIdx+0].barycenter = Eigen::Vector3f(1.0, 0.0, 0.0);
            vertices[fIdx+1].barycenter = Eigen::Vector3f(0.0, 1.0, 0.0);
            vertices[fIdx+2].barycenter = Eigen::Vector3f(0.0, 0.0, 1.0);
        }
    }
}

void GLMesh::setup(const std::vector<Eigen::Vector3f>& colors)
{
    fillBuffers(colors);
    
    // bind vao
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // bind vbo
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, (GLsizei)vertices.size() * sizeof(GLVertex), &vertices[0], GL_STATIC_DRAW);

    // set vertex attribute pointers for vbo data
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (GLvoid*)0);
    
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (GLvoid*)offsetof(GLVertex, normal));
    
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (GLvoid*)offsetof(GLVertex, barycenter));
    
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (GLvoid*)offsetof(GLVertex, color));
    
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 2, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (GLvoid*)offsetof(GLVertex, uv));
    
    // unbind vao
    glBindVertexArray(0);
}

void GLMesh::update(const std::vector<Eigen::Vector3f>& colors)
{
    fillBuffers(colors);
    
    // bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, (GLsizei)vertices.size() * sizeof(GLVertex), &vertices[0]);
}

void GLMesh::draw(Shader& shader) const
{
    shader.use();    
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, (GLsizei)vertices.size());
    glBindVertexArray(0);
}

void GLMesh::reset()
{
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
}

GLPoint::GLPoint(const Eigen::Vector3f& point0):
point(point0),
vao(0),
vbo(0)
{
    
}

void GLPoint::setup()
{
    // bind vao
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    
    // bind vbo
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Eigen::Vector3f), point.data(), GL_STATIC_DRAW);
    
    // set vertex attribute pointers for vbo data
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (GLvoid*)0);
    
    // unbind vao
    glBindVertexArray(0);
}

void GLPoint::draw(Shader& shader) const
{
    shader.use();
    glPointSize(4.0);
    glBindVertexArray(vao);
    glDrawArrays(GL_POINTS, 0, 1);
    glBindVertexArray(0);
}

void GLPoint::reset()
{
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
}
