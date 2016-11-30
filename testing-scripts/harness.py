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

        # Run each method
        for method in methods:

            ### Handle directories for this dataset and method
            datasetPath = os.path.join(dataRoot, dataset)

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




    #### Generate plots ####




if __name__ == "__main__":
    main()
