# Basic application to load a mesh from file and solve a Poisson problem on it

# Python imports
import sys, random, time
from os.path import basename

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix, csc_matrix
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
        raise Exception("Syntax: nEigs path/to/mesh.obj path/to/outfile.txt")

    # Read in the mesh
    mesh = cy.HalfEdgeMesh(cy.readMesh(filename))


    ### Compute eigenvalues
    print("=== Computing the smallest " + str(nEigs) + " eigenvectors of the Laplacian ===")


    ## Build the laplacian matrix
    print("Setting up matrix...")
    t0 = time.time()

    ind = mesh.enumerateVertices()
    N = len(mesh.verts)
    A = lil_matrix((N,N))
    M = lil_matrix((N,N))

    # Iterate over the interior vertices, adding terms for the neighbors of each
    EPS = 1e-6
    for vert in mesh.verts:

        area = vert.dualArea

        i = ind[vert]
        for (neighEdge, neighVert) in vert.adjacentEdgeVertexPairs():
            w = neighEdge.cotanLaplace

            if w < 0: w = 0.0 # ensure diagonal dominance

            j = ind[neighVert]
            A[i,j] -= w
            A[i,i] += w

        A[i,i] += EPS

        M[i,i] = area

    A = csc_matrix(A)
    M = csc_matrix(M)

    # Announce and timing
    tSetup = time.time() - t0
    print("  ...setup complete.")
    print("Calling ARPACK to solve...") 
    t0 = time.time()

    ## Solve the problem using ARPACK
    # eigenvalues, eigenvectors = scipy.sparse.linalg.eigsh(A, nEigs, M=M, which="SM")
    eigenvalues, eigenvectors = scipy.sparse.linalg.eigsh(A, nEigs, M=M, sigma=0.0) # ooh, ahhh shift-invert mode
 
    # Measure error
    # errorMat = (A.dot(eigenvectors) - M.dot(eigenvectors) * eigenvalues)
    # print(errorMat.shape)
    # error = errorMat.max()
    # print("Max error = " + str(error))

    # Normalize eigenvectors such that the first element is positive and the norm is 1
    norms = scipy.linalg.norm(eigenvectors, axis=0)
    signs = np.sign(eigenvectors[0,:])
    eigenvectors = eigenvectors * signs / norms

    # Announce and timing
    tSolve = time.time() - t0
    print("  ...solve complete.")
    print("  setup time: "+str(tSetup)+" solve time: " + str(tSolve))

    # print(eigenvalues)
    
    # Test
    # for i in range(nEigs):
        # print("\nTesting i = " + str(i))
        # testEvec = eigenvectors[:, i]
        # testEval = eigenvalues[i]
        # # diff = A.dot(testEvec) - testEval * testEvec
        # diff = A.dot(testEvec) / (testEval * testEvec)
        # print(diff)
    

    # Open output file and write result
    with open(outFilename, 'w') as outFile:

        # Write the eigenvectors
        for iVert in range(len(mesh.verts)):

            for iEig in range(nEigs):

                outFile.write(str(eigenvectors[iVert, iEig]) + " ")

            outFile.write("\n")

if __name__ == "__main__": main()
