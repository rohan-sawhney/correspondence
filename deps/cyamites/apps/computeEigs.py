# Basic application to load a mesh from file and solve a Poisson problem on it

# Python imports
import sys, random, time
from os.path import basename

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
import scipy.sparse.linalg

# Cyamites imports
import cyamites as cy

def main():


    # Get the path for the mesh to load, either from the program argument if
    # one was given, or a dialog otherwise
    if(len(sys.argv) > 3):
        nEigs = int(sys.argv[1])
        filename = sys.argv[2]
        outFilename = sys.argv[3]
    else:
        # TODO Implement a file dialog here
        raise Exception("Syntax: nEigs path/to/mesh.obj")

    # Read in the mesh
    mesh = cy.HalfEdgeMesh(cy.readMesh(filename))


    ### Compute eigenvalues
    print("=== Computing the smallest " + str(nEigs) + " eigenvalues of the Laplacian ===")


    ## Build the laplacian matrix
    print("Setting up matrix...")
    t0 = time.time()

    ind = mesh.enumerateVertices()
    N = len(mesh.verts)
    A = lil_matrix((N,N))

    # Iterate over the interior vertices, adding terms for the neighbors of each
    for vert in mesh.verts:

        i = ind[vert]
        for (neighEdge, neighVert) in vert.adjacentEdgeVertexPairs():
            w = neighEdge.cotanLaplace
            j = ind[neighVert]
            A[i,j] -= w
            A[i,i] += w


    # Announce and timing
    tSetup = time.time() - t0
    print("  ...setup complete.")
    print("Calling ARPACK to solve...") 
    t0 = time.time()

    ## Solve the problem using ARPACK
    eigenvalues, eigenvectors = scipy.sparse.linalg.eigsh(A, nEigs, which="SM")
    
    
    # Announce and timing
    tSolve = time.time() - t0
    print("  ...solve complete.")
    print("  setup time: "+str(tSetup)+" solve time: " + str(tSolve))


    # Open output file and write result
    with open(outFilename, 'w') as outFile:

        outFile.write(str(N) + "\n")
        outFile.write(str(nEigs) + "\n")

        # Write the eigenvalues
        for x in eigenvalues:
            outFile.write(str(x) + "\n")

        # Write the eigenvectors
        for vert in mesh.verts:

            for iEig in range(nEigs):

                outFile.write(str(eigenvectors[ind[vert], iEig]) + " ")

            outFile.write("\n")

if __name__ == "__main__": main()
