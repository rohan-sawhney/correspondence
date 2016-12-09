# General purpose harness to run and evaluate algorithms for shape descriptors

import os, argparse, subprocess, time
from timeit import default_timer as timer

### Wrappers to evaluate each of the algorithms.
# Should return a wall clock time in seconds in addition to writing output files

def eval_strawman(inputFile, outputFile):

    methodLocation = "../methods/strawman/strawman-features.py"
  
    # Start timing 
    startTime = timer()
    
    # Start the process
    runner = subprocess.Popen(["python", methodLocation, inputFile, outputFile])
    runner.wait()
    
    # Finish timing 
    elapsedTime = timer() - startTime

    return elapsedTime


def eval_HKS():
    pass


### Helpers
def prettyPrintTime(elapsed):

    hours, rem = divmod(elapsed, 3600)
    minutes, seconds = divmod(rem, 60)

    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds)



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

    return gList

### Primary method

def main():

    # ==== Argument parsing ====
    parser = argparse.ArgumentParser(description="General harness to run and evaluate shape descriptor algorithms")

    parser.add_argument("root", help="Root path containing data and results")
    parser.add_argument("-f", "--force", help="Overwrite existing results", action="store_true")

    args = parser.parse_args()

    # File location things
    root = args.root;
    dataRoot = os.path.join(root, "data/")
    featuresRoot = os.path.join(root, "features/")

    # Build dataset lists
    # datasets = ["TOSCA","FAUST"]
    datasets = ["TOSCA"]

    # Build method lists
    methods = ["strawman"]
    methodFunctions = {"strawman" : eval_strawman,
                       "HKS" : eval_HKS}

    # Three stages:
    # 1) Run and benchmark algorithms to generate features on the datasets
    # 2) Evaluate the goodness of each algorithm using the ground truth data
    # 3) Generate all kinds of plots


    #### Run algorithms ####
    for dataset in datasets:

        # Handle directories for this dataset and method
        datasetPath = os.path.join(dataRoot, dataset)

        # Make sure that there is at least one obj file in the dataset
        haveOBJ = False
        for f in os.listdir(datasetPath):
            if f.endswith(".obj"):
                haveOBJ = True
        if not haveOBJ:
            print("ERROR: Found 0 obj files for dataset at " + datasetPath)
            exit()


        # Run each method
        for method in methods:

            # Ensure the dataset directory exists
            if(not os.path.exists(datasetPath)):
                print("ERROR: Could not find dataset " + datasetPath)
                exit()
            
            # Create output location if needed
            outDir = os.path.join(featuresRoot, dataset, method)
            if(not os.path.exists(outDir)):
                os.makedirs(outDir)
            
            # Open a timings file, abort if it is already present
            timingsFilename = os.path.join(featuresRoot, dataset, method, "timings.txt")
            if(os.path.exists(timingsFilename)):
                if(args.force):
                    os.remove(timingsFilename)
                else:
                    print("Skipping method: " + method + " on dataset " + dataset + ". Timings file already present")
                    continue
            timingsFile = open(timingsFilename, 'w')


            ### Process each obj file for the dataset
            for f in os.listdir(datasetPath):

                fullname = os.path.join(datasetPath, f)

                if not os.path.isfile(fullname):
                    continue
                if not f.endswith(".obj"):
                    continue

                basename = os.path.basename(f)
                name, _ = os.path.splitext(basename)
                fullnameNoExt, _ = os.path.splitext(fullname)

                inFile = fullname

                print("\nProcessing file " + str(name))
                print("  infile " + str(inFile))

                print("  Using method " + method)
           
                outFile = os.path.join(outDir, name + ".features")
                print("    outfile " + str(outFile))
           
                elapsedTime = methodFunctions[method](inFile, outFile)
                print("    elapsed time: " + prettyPrintTime(elapsedTime))

                timingsFile.write(name + "," + str(elapsedTime) + "\n")




    #### Evaluate accuracy ####

    # TODO This code assumes that the ground truth correspondence has
    # the implicit form of equally-numbered verticies matching up. This is
    # just one possible way for a correspondence to be represented, and should
    # be generalized.

    for dataset in datasets:

        # Handle directories for this dataset and method
        datasetPath = os.path.join(dataRoot, dataset)

        # Parse the file that tells us how the dataset is partitioned
        # in to corresponding groups
        groupFilename = os.path.join(datasetPath, "partitions.txt")
        groupList = parseGroupFile(groupFilename)

        # Evaluate each method
        for method in methods:

            # TODO generalize this for datasets that have



    #### Generate plots ####




if __name__ == "__main__":
    main()
