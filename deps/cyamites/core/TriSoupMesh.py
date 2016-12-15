import numpy as np

# Cyamites imports
from Utilities import normalize, norm

# A "triangle soup" mesh, with an array of vertices and list of face indices
# Generally, this class exists only for interfacing with external components, such as IO and visualization. Cyamites encourages
# all mesh computation to be performed on a halfedge mesh.
class TriSoupMesh(object):

    ### Construct a new mesh froms a vert array and a tri list
    def __init__(self, verts, tris, faceAttr=dict(), vertAttr=dict()):

        # A Nx3 list of vertex positions
        self.verts = np.array(verts)
        self.tris = np.array(tris, dtype=np.uint32)

        # A Mx3 list of tuples, where each gives an (i,j,k) CCW triangular face
        self.nVerts = len(verts)
        self.nTris = len(tris)

        # Dictionaries to hold various data on the faces and vertices
        # (as indexed lists): vertAttr['value'] = [val1, val2, val3...]
        self.faceAttr = faceAttr
        self.vertAttr = vertAttr

        # Numpy-ify all the attr data
        # TODO really need some kind of more general solution here... this
        # will fail on lots of things
        for key in self.faceAttr:
            self.faceAttr[key] = np.array(self.faceAttr[key])
        for key in self.vertAttr:
            self.vertAttr[key] = np.array(self.vertAttr[key])


    # Compute face and vertex normals, using numpy operations for efficiency
    # This will generally be much faster than computing normals on the halfedge
    # mesh.
    # 'useExisting' will cause this method to exit immediately if there is already
    #               a 'normal' attribute defined on the vertices
    def computeNormals(self, useExisting=False):

        if useExisting and 'normal' in self.vertAttr:
            #print("Skipping normal computation, normals are already defined")
            return

        # Expand out an indexed view of our vertices
        faceVerts = self.verts[self.tris]

        # Compute face normals as a cross product of the face vertices
        # TODO why the double colon...?
        faceNormals = np.cross(faceVerts[::,1 ] - faceVerts[::,0],
                               faceVerts[::,2 ] - faceVerts[::,0])

        # Area of each face is used for a weighted average below, so save that
        faceAreas = 0.5 * norm(faceNormals, axis=1).reshape((self.nTris, 1))
        faceNormals = normalize(faceNormals)

        self.faceAttr['normal'] = faceNormals


        # Now compute vertex normals as a area-weighted average of the face normals
        vertNormals = np.zeros( self.verts.shape, dtype=self.verts.dtype)
        vertNormals[ self.tris[:,0] ] += faceNormals * faceAreas
        vertNormals[ self.tris[:,1] ] += faceNormals * faceAreas
        vertNormals[ self.tris[:,2] ] += faceNormals * faceAreas
        normalize(vertNormals)

        self.vertAttr['normal'] = vertNormals

    def computeCurvature(self):
        pass
