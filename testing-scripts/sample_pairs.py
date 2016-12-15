# Read a partitions file and emit a pairs file by subsampling pairs

import sys, random

def parseGroupFile(filename):

    print("Parsing group file " + filename)

    groupLists = {}

    with open(filename) as groupFile:
        while True:

            try:
                line = next(groupFile).strip()
            except:
                break

            if(line[0] == '#'): continue

            # Must start with "GROUP:"
            line = line[6:]

            # Parse out name and count
            items = line.split(",")
            gName = items[0]
            gCount = int(items[1])

            # Accumulate the entries in the group
            gList = []
            for i in range(gCount):
                line = next(groupFile)
                gList.append(line.strip())
                
            groupLists[gName] = gList

    print("Found " + str(len(groupLists)) + " groups")

    return groupLists


def main():

    if(len(sys.argv) == 4):
        nPairsPer = int(sys.argv[1])
        partitionFile = sys.argv[2]
        pairFile = sys.argv[3]
    else:
        # TODO Implement a file dialog here
        raise Exception("Syntax: nPairsPer inPartitionsFile outPairsFile")



    # Read in the group file
    groups = parseGroupFile(partitionFile)


    # Sample pairs
    pairs = []
    for groupName in groups:     

        group = groups[groupName]

        for itemA in group:

            # all of the elements except itemA
            restOfGroup = list(group)
            restOfGroup.remove(itemA)

            # how many to sample
            sampleCount = min([nPairsPer, len(restOfGroup)])
            
            otherItems = random.sample(restOfGroup, sampleCount)
   
            for itemB in otherItems:

                pairs.append((itemA, itemB))




    # Write the pairs to a file
    with open(pairFile, 'w') as outFile:
        
        for pair in pairs:

            outFile.write(pair[0] + "," + pair[1] + "\n")

if __name__ == "__main__":
    main()
