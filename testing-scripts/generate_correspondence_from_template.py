# Basic application to load a mesh from file and solve a Poisson problem on it

# Python imports
import sys, random, time, os
from os.path import basename

import numpy as np
import scipy.spatial

# Cyamites imports
import cyamites as cy




def main():

    if(len(sys.argv) == 5):
        pairsFile = sys.argv[1]
        meshDirectory = sys.argv[2]
        templateDirectory = sys.argv[3]
        outDirectory = sys.argv[4]
    else:
        # TODO Implement a file dialog here
        raise Exception("Syntax: pairsFile meshDir templateDir")


    # Process each pair
    oldMeshAName = ""
    for line in open(pairsFile):

        if(line[0] == "#"): continue
       
        # Parse the pair
        items = line.strip().split(",")
        meshAName = items[0]
        meshBName = items[1]
        print("\nRunning on pair: " + str(meshAName) + " -- " + str(meshBName))

       
        # Build relevant filenames
        meshAFilename = os.path.join(meshDirectory, meshAName + ".obj")
        meshBFilename = os.path.join(meshDirectory, meshBName + ".obj")
        templateAFilename = os.path.join(templateDirectory, meshAName.replace("scan","reg") + ".ply")
        templateBFilename = os.path.join(templateDirectory, meshBName.replace("scan","reg") + ".ply")
        outFilename = os.path.join(outDirectory, meshAName + "-" + meshBName + ".groundtruth")
    
        # Read in the meshes
        if meshAName != oldMeshAName: # hack to do half as many mesh parses
            meshA = cy.readMesh(meshAFilename)
            templateA = cy.readMesh(templateAFilename)
        meshB = cy.readMesh(meshBFilename)
        templateB = cy.readMesh(templateBFilename)

        # Do work
        runOnPair(meshA, templateA, meshB, templateB, outFilename)

        oldMeshAName = meshAName

def runOnPair(meshA, templateA, meshB, templateB, outFilename):


    # Build nearest-neighbor structures for both meshes
    aTree = scipy.spatial.KDTree(meshA.verts) 
    bTree = scipy.spatial.KDTree(meshB.verts) 
  

    # For each point on each of the templates, find the nearest point on the other shape
    def mapPoints(sourceMesh, targetKDtree):

        targetResults = []

        for pos in sourceMesh.verts:
   
            _, ind = targetKDtree.query(pos)
            targetResults.append(ind)

        return targetResults

    aNearests = mapPoints(templateA, aTree)
    bNearests = mapPoints(templateB, bTree)


    # Open output file and write out the pairs of matching points
    with open(outFilename, 'w') as outFile:

        # Write the eigenvectors
        for aInd,bInd in zip(aNearests, bNearests):
            outFile.write(str(aInd) + "," + str(bInd) + "\n")

            # print(str(meshA.verts[aInd])  + " -- " + str(meshB.verts[bInd]))

if __name__ == "__main__": main()
