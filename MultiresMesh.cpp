#include "MultiresMesh.h"
#include "Mesh.h"
#include "MeshIO.h"
#define _USE_MATH_DEFINES
#include <openmesh/Tools/Decimater/DecimaterT.hh>
#include <openmesh/Tools/Decimater/ModQuadricT.hh>
#undef _USE_MATH_DEFINES

typedef OpenMesh::Decimater::DecimaterT<OMesh> ODecimater;
typedef OpenMesh::Decimater::ModQuadricT<OMesh>::Handle OModQuadric;

MultiresMesh::MultiresMesh(Mesh *mesh0, double d0, double C0):
d(d0),
C(C0)
{
    lods.push_back(mesh0);
    buildOMesh(lods[0]);
}

MultiresMesh::~MultiresMesh()
{
    for (int i = 1; i < (int)lods.size(); i++) {
        delete lods[i];
    }
}

void MultiresMesh::buildOMesh(Mesh *mesh)
{
    // set vertices
    std::vector<OMesh::VertexHandle> vertexHandles(mesh->vertices.size());
    for (VertexCIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        vertexHandles[v->index] = oMesh.add_vertex(OMesh::Point(v->position.x(),
                                                                v->position.y(),
                                                                v->position.z()));
    }
    
    // set faces
    for (FaceCIter f = mesh->faces.begin(); f != mesh->faces.end(); f++) {
        OMesh::FaceHandle faceHandle = oMesh.add_face(vertexHandles[f->he->vertex->index],
                                                      vertexHandles[f->he->next->vertex->index],
                                                      vertexHandles[f->he->next->next->vertex->index]);
        if (!faceHandle.is_valid()) {
            std::cout << "Invalid face handle: " << f->index << std::endl;
        }
    }
}

void MultiresMesh::buildMesh(Mesh *mesh)
{
    MeshData data;
    
    // set vertices
    for (OMesh::ConstVertexIter v_it = oMesh.vertices_begin(); v_it != oMesh.vertices_end(); ++v_it) {
        const OMesh::Point& point(oMesh.point(*v_it));
        data.positions.push_back(Eigen::Vector3d(point[0], point[1], point[2]));
    }
    
    // set faces
    for (OMesh::ConstFaceIter f_it = oMesh.faces_begin(); f_it != oMesh.faces_end(); ++f_it) {
        // set indices
        std::vector<Index> faceIndices;
        OMesh::ConstFaceVertexIter fv_it = oMesh.fv_iter(*f_it);
        faceIndices.push_back(Index(fv_it->idx(), -1, -1)); ++fv_it;
        faceIndices.push_back(Index(fv_it->idx(), -1, -1)); ++fv_it;
        faceIndices.push_back(Index(fv_it->idx(), -1, -1));
        
        data.indices.push_back(faceIndices);
    }
    
    MeshIO::buildMesh(data, *mesh);
}

Mesh* MultiresMesh::decimate(int v)
{
    // add modules to decimator
    ODecimater oDecimater(oMesh);
    
    OModQuadric oModQuadric;
    oDecimater.add(oModQuadric);
    
    if (oDecimater.initialize()) {
        // decimate
        oDecimater.decimate_to(v);
        
        // clean up
        oMesh.garbage_collection();
        
        // build mesh
        Mesh *mesh = new Mesh();
        buildMesh(mesh);
        
        return mesh;
    }
    
    std::cout << "Unable to initialize decimater" << std::endl;
    return NULL;
}

void MultiresMesh::build()
{
    int n = (int)lods[0]->vertices.size();
    while (n > C) {
        n = n / d;
        
        // decimate
        lods.push_back(decimate(n));
        
        // TODO: project
    }
}
