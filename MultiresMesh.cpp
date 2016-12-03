#include "MultiresMesh.h"
#include "Mesh.h"
#include "MeshIO.h"
#include "Bvh.h"
#include <future>
#include <functional>
#define _USE_MATH_DEFINES
#include <openmesh/Tools/Decimater/DecimaterT.hh>
#include <openmesh/Tools/Decimater/ModQuadricT.hh>
#include <OpenMesh/Tools/Decimater/ModNormalFlippingT.hh>
#undef _USE_MATH_DEFINES

typedef OpenMesh::Decimater::DecimaterT<OMesh> ODecimater;
typedef OpenMesh::Decimater::ModQuadricT<OMesh>::Handle OModQuadric;
typedef OpenMesh::Decimater::ModNormalFlippingT<OMesh>::Handle OModNormalFlipping;

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
    oMesh.request_face_normals();
    
    // set vertices
    std::vector<OMesh::VertexHandle> vertexHandles(mesh->vertices.size());
    for (VertexCIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        vertexHandles[v->index] = oMesh.add_vertex(OMesh::Point(v->position.x(),
                                                                v->position.y(),
                                                                v->position.z()));
    }
    
    // set faces
    for (FaceCIter f = mesh->faces.begin(); f != mesh->faces.end(); f++) {
        if (!f->isBoundary()) {
            OMesh::FaceHandle faceHandle = oMesh.add_face(vertexHandles[f->he->vertex->index],
                                                          vertexHandles[f->he->next->vertex->index],
                                                          vertexHandles[f->he->next->next->vertex->index]);
            if (!faceHandle.is_valid()) {
                std::cout << "Invalid face handle: " << f->index << std::endl;
            }
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
        OMesh::ConstFaceVertexCCWIter fv_it = oMesh.fv_ccwiter(*f_it);
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
    
    OModNormalFlipping oModNormalFlipping;
    oDecimater.add(oModNormalFlipping);
    oDecimater.module(oModNormalFlipping).set_max_normal_deviation(30.0);
    
    if (oDecimater.initialize()) {
        // update face normals
        oMesh.update_face_normals();
        
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

void MultiresMesh::project(Mesh *mesh1, Mesh *mesh2)
{
    // build bvh
    Bvh bvh(mesh2);
    bvh.build();
    
    // project vertices
    std::vector<std::future<int>> futures; futures.reserve(mesh1->vertices.size());
    for (VertexIter v = mesh1->vertices.begin(); v != mesh1->vertices.end(); v++) {
        v->projection.d = INFINITY;
        futures.push_back(std::async(std::launch::async, &Bvh::getIntersection,
                                     &bvh, std::ref(v->projection.d),
                                     std::ref(v->projection.p), std::cref(v->position)));
    }
    
    for (VertexIter v = mesh1->vertices.begin(); v != mesh1->vertices.end(); v++) {
        v->projection.fIdx = futures[v->index].get();
    }
}

void MultiresMesh::build()
{
    int i = 0;
    int n = (int)lods[0]->vertices.size();
    while (n > C) {
        n = n / d;
        
        // decimate
        Mesh *mesh = decimate(n);
        if (mesh) lods.push_back(mesh);
        else break;
        
        // project
        project(lods[i], lods[i+1]);

        i++;
    }
}

int MultiresMesh::numLods()
{
    return (int)lods.size();
}

Mesh* MultiresMesh::lod(int l)
{
    if (l >= 0 && l < (int)lods.size()) {
        return lods[l];
    }
    
    return NULL;
}
