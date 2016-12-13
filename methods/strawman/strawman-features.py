# A simple dummy program which uses the x,y,z coordinates of each vertex as features

import argparse

def readMesh_OBJ(filename):

    verts = []
    tris = []

    # Process the file line by line
    lineInd = -1
    for line in open(filename).readlines():
        lineInd += 1
        if line[0] == '#':
            continue

        items = line.strip().split(' ')

        # Process vertex
        if items[0] == 'v':
            verts.append([float(s) for s in items[1:]])

        # Process tex-coord
        elif items[0] == 'vt':
            pass

        # Process normal vector
        elif items[0] == 'vn' : #normal vector
            pass

        # Process face indices
        elif items[0] == 'f' :
            face = items[1:]
            if len(face) != 3 :
                print("Line " + str(lineInd) + " is not a triangular face: " + line)
                continue

            # Take index before slash, which is the vertex index (if slash is given).
            # Conversion from 1-based to 0-based indexing happens at end
            tris.append([int(s.split("/")[0]) for s in face])

    # Convert to numpy arrays
    # verts = np.array(verts)
    # tris = np.array(tris)

    # Convert to 0-based indexing
    # tris = tris - 1

    return verts, tris

def main():


    # ==== Argument parsing ====
    parser = argparse.ArgumentParser(description="Compute coordinate features")

    parser.add_argument("input", help="Input file")
    parser.add_argument("output", help="Output file")

    args = parser.parse_args()

    inputFile = args.input
    outputFile = args.output

    # Parse input
    verts, tris = readMesh_OBJ(inputFile)


    # Write output
    with open(outputFile, 'w') as outFile:

        for v in verts:

            outFile.write(str(v[0]) + " " + str(v[1]) + " " + str(v[2]) + "\n")



if __name__ == "__main__":
    main()
