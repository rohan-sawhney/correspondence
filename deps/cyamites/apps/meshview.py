# Basic application to load a mesh from file and view it in window
# using openGL

# Python imports
import sys
from os.path import basename
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

    # Toss up a viewer window
    winName = 'Cyamites meshview -- ' + basename(filename)
    meshDisplay = cy.MeshDisplay(windowTitle=winName)
    meshDisplay.setMesh(mesh)
    meshDisplay.setVectors("normal", vectorDefinedAt='face')
    meshDisplay.startMainLoop()


if __name__ == "__main__": main()
