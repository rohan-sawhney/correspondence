# Test displaying a uniform vector field on a mesh

# Python imports
import sys, time
from os.path import basename
from math import *
import numpy as np

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

    # Uniform direction
    unifDir = cy.normalize(np.array([0.0,0.0,1.0]))

    # For each vertex, set a reference edge/direction and project the uniform
    # direction in to the normal plane
    print("Computing test directions...")

    startTime = time.time()
    mesh.assignReferenceDirections()
    for vert in mesh.verts:

        # Project
        uNorm = np.dot(vert.normal, unifDir)*vert.normal
        uInPlaneDir = cy.normalize(unifDir - uNorm)
        xVal = np.dot(vert.referenceDirectionR3, uInPlaneDir)
        yVal = np.dot(np.cross(vert.normal, vert.referenceDirectionR3), uInPlaneDir)
        vert.theta = atan2(yVal, xVal)

    elapsedTime = time.time() - startTime
    print("Computing directions took " + str(elapsedTime))

    # Toss up a viewer window
    winName = 'Cyamites meshview -- ' + basename(filename)
    meshDisplay = cy.MeshDisplay(windowTitle=winName)
    meshDisplay.setVectors('theta', nSym=1, isTangentVector=True, isUnit=True)
    meshDisplay.setMesh(mesh)
    meshDisplay.startMainLoop()


if __name__ == "__main__": main()
