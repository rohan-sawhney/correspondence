# Basic application to load a mesh from file and solve a Poisson problem on it

# Python imports
import sys, random
from os.path import basename

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve

# PyAMG algebraic multigrid
from pyamg import *

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
    mesh = cy.HalfEdgeMesh(cy.readMesh(filename), staticGeometry=False)

    # Check for degenerate faces, because those will break curvature flow
    mesh.checkDegenerateFaces()

    # Create the solver object
    curvatureFlower = cy.MeanCurvatureFlower(mesh, hStep = 1.0)


    # Toss up a viewer window
    winName = 'Cyamites Mean Curvature Flow -- ' + basename(filename)
    meshDisplay = cy.MeshDisplay(windowTitle=winName)
    meshDisplay.setMesh(mesh)


    # Add forward and backward position stacks to each vertex (think browser history)
    for v in mesh.verts:
        v.backstack = []
        v.forwardstack = []
    stackSizes = [0,0] # This is a list partly to hack around weird python scoping

    # Define a callback in the viewer which flows forward, either computing a new
    # step of flow or popping from the forward stack if we've already saved the next
    # position.
    def flowForward():

        if stackSizes[1] == 0:
            for v in mesh.verts:
                v.backstack.append(v.pos.copy())
            curvatureFlower.flowStep()
            stackSizes[0] += 1
        else:
            stackSizes[1] -= 1
            stackSizes[0] += 1
            for v in mesh.verts:
                v.backstack.append(v.pos)
                v.pos = v.forwardstack.pop()

        meshDisplay.readNewMeshValues()

    # Flows backwards if possible, does nothing otherwise
    def flowBackward():

        if stackSizes[0] > 0:
            stackSizes[0] -= 1
            stackSizes[1] += 1
            for v in mesh.verts:
                v.forwardstack.append(v.pos.copy())
                v.pos = v.backstack.pop()

        meshDisplay.readNewMeshValues()

    meshDisplay.registerKeyCallback('x', flowForward, "Perform one step of mean curvature flow")
    meshDisplay.registerKeyCallback('z', flowBackward, "Undo one step of mean curvature flow")

    meshDisplay.startMainLoop()


if __name__ == "__main__": main()
