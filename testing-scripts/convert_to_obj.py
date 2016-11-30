# Simple script for converting between some geometry file formats.

import argparse
import os

def parseVertTri(vertFilename, triFilename):

    # Read verts
    verts = []
    with open(vertFilename) as vertFile:
        for line in vertFile:
            items = line.strip().split(" ")
            v = (float(items[0]), float(items[1]), float(items[2]))
            verts.append(v)


    # Read tris
    tris = []
    with open(triFilename) as triFile:
        for line in triFile:
            items = line.strip().split(" ")
            f = (int(items[0]), int(items[1]), int(items[2]))
            tris.append(f)


    return verts, tris

def writeOBJ(verts, tris, filename, force=False):

    # Check if the file exists
    if(not force and os.path.exists(filename)):
        print("Not writing output; file exists at " + filename + " .  Use --force to override")
        return

    print("Writing output to " + filename)

    with open(filename, 'w') as outFile:

        for v in verts:
            outFile.write("v " + str(v[0]) + " " + str(v[1]) + " " + str(v[2]) + "\n")
        
        for f in tris:
            outFile.write("f " + str(f[0]) + " " + str(f[1]) + " " + str(f[2]) + "\n")


def main():

    # ==== Argument parsing ====
    parser = argparse.ArgumentParser(description="Convert mesh filetypes")

    parser.add_argument("input", nargs='*', help="Input file")

    parser.add_argument("-f", "--force", help="Overwrite existing output file", action="store_true")
    args = parser.parse_args()

    inputPathList = args.input
    overwrite = args.force 

    for inputPath in inputPathList: 

        inputDir, inputName = os.path.split(inputPath)
        inputPrefix, inputExt = os.path.splitext(inputPath)
        outputName = inputPrefix + ".obj" 


        print("Processing file " + inputPath)


        # Dispatch to the proper version
        if(inputExt == ".vert"):
            vertFilename = inputPath
            triFilename = inputPrefix + ".tri"
            verts, tris = parseVertTri(vertFilename, triFilename) 

        else:
            print("Extension " + inputExt + " not recognized. Exiting")
            exit()


        writeOBJ(verts, tris, outputName, overwrite)


if __name__ == "__main__":
    main()
