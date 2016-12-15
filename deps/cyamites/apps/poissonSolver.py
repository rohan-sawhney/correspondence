# Basic application to load a mesh from file and solve a Poisson problem on it

# Python imports
import sys, random
from os.path import basename

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve

# Cyamites imports
import cyamites as cy

def main():


    # Get the path for the mesh to load, either from the program argument if
    # one was given, or a dialog otherwise
    if(len(sys.argv) > 1):
        filename = sys.argv[1]
    else:
        # TODO Implement a file dialog here
        raise Exception("No file name specified")

    # Read in the mesh
    mesh = cy.HalfEdgeMesh(cy.readMesh(filename))


    ### Solve a poisson problem (heat equation)
    print("Building Poisson problem")

    # Pick some vertices as boundary values and assign them a random value
    nBoundaryVerts = 15
    boundaryVerts = set(random.sample(mesh.verts, nBoundaryVerts))
    boundaryValue = dict()
    for v in boundaryVerts:
        boundaryValue[v] = random.random()

    # Solve the Poisson problem
    sol = cy.solvePoisson(mesh, boundaryValue, method='pyamg-RS')
    mesh.applyVertexValue(sol, "heatVal")


    # Toss up a viewer window
    print("Viewing solution")
    winName = 'Cyamites Poisson Solver -- ' + basename(filename)
    meshDisplay = cy.MeshDisplay(windowTitle=winName)
    meshDisplay.setMesh(mesh)
    meshDisplay.setShapeColorFromScalar("heatVal", cmapName="seismic", vMinMax=[0.0,1.0])
    meshDisplay.startMainLoop()

if __name__ == "__main__": main()
