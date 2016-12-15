# Basic application to load a mesh from file and solve a Poisson problem on it

# Python imports
import sys, random, time
from os.path import basename

import numpy as np

# Cyamites imports
import cyamites as cy

def main():


    # Get the path for the mesh to load, either from the program argument if
    # one was given, or a dialog otherwise
    if(len(sys.argv) > 7):
        meshAFilename = sys.argv[2]
        templateAFilename = sys.argv[3]
        meshBFilename = sys.argv[4]
        templateBFilename = sys.argv[5]
        outFilename = sys.argv[6]
    else:
        # TODO Implement a file dialog here
        raise Exception("Syntax: nEigs path/to/mesh.obj path/to/outfile.txt")

    # Read in the meshes
    meshA = cy.HalfEdgeMesh(cy.readMesh(meshAFilename))
    templateA = cy.HalfEdgeMesh(cy.readMesh(templateAFilename))
    meshB = cy.HalfEdgeMesh(cy.readMesh(meshBFilename))
    templateB = cy.HalfEdgeMesh(cy.readMesh(templateBFilename))

    # Create numpy arrays of positions for each mesh
    def vertPosToArray(mesh):

        pList = []
        for v in mesh.verts:
            pList.append(v.pos)

        return np.array(pList) 

    meshAVerts = vertPosToArray(meshA)
    templateAVerts = vertPosToArray(templateA)
    meshBVerts = vertPosToArray(meshB)
    templateBVerts = vertPosToArray(templateB)


    # Open output file and write result
    with open(outFilename, 'w') as outFile:

        # Write the eigenvectors
        for iVert in range(len(mesh.verts)):

            for iEig in range(nEigs):

                outFile.write(str(eigenvectors[iVert, iEig]) + " ")

            outFile.write("\n")

if __name__ == "__main__": main()
