# General purpose harness to run and evaluate algorithms for shape descriptors

import os, argparse, subprocess, time, random
from timeit import default_timer as timer

import numpy as np
import scipy
import scipy.interpolate
import matplotlib.pyplot as plt
import seaborn as sns

# Make the args globally available
global args

# Number of precision increments at which to compute precision-recall curve value
nPRIncrements = 500

### Wrappers to compute the result each of the algorithms.
# Should return a wall clock time in seconds in addition to writing output files

def compute_spectrum(n, evalsAndEvecs, inputFile, outputFile):

    methodLocation = "../methods/laplacian-eigenvec/computeEVecs.py"
    
    # Start timing
    startTime = timer()
    
    # Start the process
    runner = subprocess.Popen(["python", methodLocation, n, evalsAndEvecs, inputFile, outputFile])
    runner.wait()
    
    # Finish timing
    elapsedTime = timer() - startTime
    
    return elapsedTime


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


def eval_laplacian(inputFile, outputFile):

    return compute_spectrum("30", "n", inputFile, outputFile)


def eval_curvature(inputFile, outputFile):

    methodLocation = "../build/correspond"
  
    # Start timing 
    startTime = timer()
    
    # Start the process
    runner = subprocess.Popen([methodLocation, "-descriptor", "3", "-obj_path", inputFile, "-output_path", outputFile])
    runner.wait()
    
    # Finish timing 
    elapsedTime = timer() - startTime

    return elapsedTime


def eval_HKS(inputFile, outputFile, eigFile):

    methodLocation = "../build/correspond"
    
    # Start timing
    startTime = timer()
    
    # Start the process
    runner = subprocess.Popen([methodLocation, "-descriptor", "0", "-obj_path", inputFile,
                               "-eig_path", eigFile, "-output_path", outputFile])
    runner.wait()
    
    # Finish timing
    elapsedTime = timer() - startTime
    
    return elapsedTime


def eval_WKS(inputFile, outputFile, eigFile):

    methodLocation = "../build/correspond"
    
    # Start timing
    startTime = timer()
    
    # Start the process
    runner = subprocess.Popen([methodLocation, "-descriptor", "2", "-obj_path", inputFile,
                               "-eig_path", eigFile, "-output_path", outputFile])
    runner.wait()
    
    # Finish timing
    elapsedTime = timer() - startTime
    
    return elapsedTime

### Evaluate the accuracy of a method

# Parses a feature file in to a numpy array
def parseFeatureFile(featureFile):

    return np.loadtxt(featureFile)


# Takes two same-shaped input arrays, each representing a mesh, where the i'th row is
# the feature vector for that vertex. It is assumed that the ground truth has the
# rows in correspondence.
# Returns a list of accuracy scores [0,1] with one value for each vertex in A. 0 is the best
# result, if its ground truth had the most similar score, and 1 is the worst if it had
# the most different score.
def evaluateAccuracy(featuresA, featuresB, subsampleRate = 1.0):

    result = []

    nRows, nCols = featuresA.shape

    for iRow in range(nRows):

        if random.random() > subsampleRate: continue
    
        # Compute distances
        iFeatures = featuresA[iRow,:]
        dists = np.linalg.norm(featuresB - iFeatures, axis=1)

        # Find the index at which the matched vertex would lie, after sorting
        sortedDists = np.argsort(dists)
        sortedInd = np.where(sortedDists == iRow)[0][0]
        relativeInd = float(sortedInd) / nRows

        result.append(relativeInd)

    return result


# Performs all appropriate pairwise calls to evaluateAccuracy and returns
# the merged result
def evaluateAllAccuracies(groupSet, featureDirectory, accuracyDirectory):


    for groupName in groupSet:
   
        print("Evaluating group: " + groupName)

        group = groupSet[groupName]

        for iGroup in range(len(group)):

            iFilename = os.path.join(featureDirectory, group[iGroup] + ".features")
            

            # Evaluate the how well we matched the corresponding vertices for each pair in the group
            for jGroup in range(len(group)):

                if iGroup == jGroup: continue
                print("Evaluating pair: " + group[iGroup] + " and " + group[jGroup])
            
                jFilename = os.path.join(featureDirectory, group[jGroup] + ".features")

                # The path to put the result at
                resultFilename = os.path.join(accuracyDirectory, group[iGroup] + "-" + group[jGroup] + ".accuracy")
                
                # Skip computation if the file exists 
                if(os.path.exists(resultFilename)):
                    if(args.force):
                        os.remove(resultFilename)
                    else:
                        print("Skipping accuracy computation for " + group[iGroup] + "---" + group[jGroup] + ", file already present")
                        continue


                # Assess the result
                # NOTE THAT WE ARE SUBSAMPLING BY A LOT RIGHT NOW
                pairwiseResult = evaluateAccuracy(parseFeatureFile(iFilename), parseFeatureFile(jFilename), 0.001)

                # Compute precision-recall samples 
                pairwisePR = computePrecisionRecall(pairwiseResult)
           
                # Write the samples to file
                np.savetxt(resultFilename, pairwisePR[1], fmt="%.8f")
                 


def parseAccuracies(groupSet, accuracyDirectory):
  
    print("Reading all accuracies from  " + accuracyDirectory)

    allAccuracies= []

    for groupName in groupSet:
   
        group = groupSet[groupName]

        for iGroup in range(len(group)):

            # Evaluate the how well we matched the corresponding vertices for each pair in the group
            for jGroup in range(len(group)):

                if iGroup == jGroup: continue

                resultFilename = os.path.join(accuracyDirectory, group[iGroup] + "-" + group[jGroup] + ".accuracy")
                acc = np.loadtxt(resultFilename)
                
                allAccuracies.append(acc)

    return allAccuracies


# Compute the value of the PR curve at nPRIncrements points between [0,1]
def computePrecisionRecall(relativeIndexVals):

    rankThreshold = np.array(sorted(relativeIndexVals))
    recallPercent = np.linspace(0, 1, num=len(relativeIndexVals))

    # Ensure we have the full range [0,1]
    # TODO ugly use of appends
    rankThreshold = np.append(np.append(np.array([0]), rankThreshold), np.array([1]))
    recallPercent = np.append(np.append(np.array([0]), recallPercent), np.array([1]))

    # Resample to regular values
    interpolator = scipy.interpolate.interp1d(rankThreshold, recallPercent)
    rankSamples = np.linspace(0, 1, num=nPRIncrements)
    recallSamples = interpolator(rankSamples)

    return (rankSamples, recallSamples)

### Plotters

def plotMethodAccuracy(methodName, accuracyLists, plotDir):

    print("Plotting accuracy for method " + methodName)

    yData = np.array(accuracyLists)
    xData = np.linspace(0, 1, num=nPRIncrements)

    sns.set(font_scale=2) 
    plt.title(methodName + " accuracy")
    plt.xlabel("Threshold")
    plt.ylabel("Fraction under threshold")

    sns.tsplot(yData, xData, err_style="unit_traces", err_palette=sns.color_palette(["#ccd8ff"]))
    # sns.tsplot(yData, xData, ci = [68,95,99.99])
    # sns.tsplot(yData, xData, err_style="boot_traces", n_boot=500)

    outname = plotDir + methodName + "-accuracy.pdf"
    plt.savefig(outname)
    plt.clf()

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

    return groupLists


### Primary method

def main():

    # ==== Argument parsing ====
    parser = argparse.ArgumentParser(description="General harness to run and evaluate shape descriptor algorithms")

    parser.add_argument("root", help="Root path containing data and results")
    parser.add_argument("-f", "--force", help="Overwrite existing results", action="store_true")

    global args
    args = parser.parse_args()

    # File location things
    root = args.root;
    dataRoot = os.path.join(root, "data/")
    featuresRoot = os.path.join(root, "features/")
    evaluateRoot = os.path.join(root, "accuracy/")
    plotRoot = os.path.join(root, "plots/")

    # Build dataset lists
    # datasets = ["TOSCA","FAUST"]
    datasets = ["TOSCA"]

    # Build method lists
    methods = ["strawman", "curvature", "laplacian-eigenvecs", "hks", "wks"]
    methodFunctions = {"strawman" : eval_strawman,
                       "curvature" : eval_curvature,
                       "laplacian-eigenvecs" : eval_laplacian,
                       "hks" : eval_HKS,
                       "wks" : eval_WKS}

    # Three stages:
    # 1) Run and benchmark algorithms to generate features on the datasets
    # 2) Evaluate the goodness of each algorithm using the ground truth data
    # 3) Generate all kinds of plots

    #### Run algorithms ####
    for dataset in datasets:

        # Handle directories for this dataset and method
        datasetPath = os.path.join(dataRoot, dataset)
        
        # Ensure the dataset directory exists
        if(not os.path.exists(datasetPath)):
            print("ERROR: Could not find dataset " + datasetPath)
            exit()

        # Make sure that there is at least one obj file in the dataset
        haveOBJ = False
        for f in os.listdir(datasetPath):
            if f.endswith(".obj"):
                haveOBJ = True
        if not haveOBJ:
            print("ERROR: Found 0 obj files for dataset at " + datasetPath)
            exit()

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
                
            # Run each method
            eigFile = ""
            eigTime = 0
            for method in methods:
                
                print("  Using method " + method)
    
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

                if (method == "hks"):
                    print("Computing Spectrum")
                    eigFile = os.path.join(outDir, name + ".eig")
                    print("    eigfile " + str(eigFile))
                    eigTime = compute_spectrum("500", "y", inFile, eigFile)

                outFile = os.path.join(outDir, name + ".features")
                print("    outfile " + str(outFile))

                if (method == "hks" or method == "wks"):
                    elapsedTime = eigTime + methodFunctions[method](inFile, outFile, eigFile)
                else:
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

            featuresDir = os.path.join(featuresRoot, dataset, method)
            evaluteDir = os.path.join(evaluateRoot, dataset, method)
            
            # Create output directory as needed
            if(not os.path.exists(evaluteDir)):
                os.makedirs(evaluteDir)

            evaluateAllAccuracies(groupList, featuresDir, evaluteDir)




    #### Generate plots ####
   
    # Make an accuracy plot for each dataset/method
    for dataset in datasets:

        # Handle directories for this dataset and method
        datasetPath = os.path.join(dataRoot, dataset)

        # Parse the file that tells us how the dataset is partitioned
        # in to corresponding groups
        groupFilename = os.path.join(datasetPath, "partitions.txt")
        groupList = parseGroupFile(groupFilename)

        methodAccuracies = {}

        for method in methods:

            evaluteDir = os.path.join(evaluateRoot, dataset, method)
   
            # Read in the evaluation results     
            evalResults = parseAccuracies(groupList, evaluteDir)

            # Make a plot of the evaluation results
            plotMethodAccuracy(method, evalResults, plotRoot)       

            # TODO ROHAN: Accumulate to make a merged plot

            # TODO ROHAN: Performance vs accuracy - area under curve


        # Make a combined plot for all methods


if __name__ == "__main__":
    main()
