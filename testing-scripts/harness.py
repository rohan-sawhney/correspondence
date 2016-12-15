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

    # NOTE: Duplicated below
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


def evaluateAccuracyWithFile(featuresA, featuresB, correspondenceFile, subsampleRate = 1.0):

    result = []

    # Read in the correspondence file
    correspondences = []
    with open(correspondenceFile) as cFile:

        for line in cFile:

            if(line[0] == "#"): continue

            items = line.strip().split(",")
            indA = int(items[0])
            indB = int(items[1])

            correspondences.append((indA,indB))



    nRows, nCols = featuresA.shape

    # NOTE: Duplicated above
    for indA, indB in correspondences:

        if random.random() > subsampleRate: continue
    
        # Compute distances
        iFeatures = featuresA[indA,:]
        dists = np.linalg.norm(featuresB - iFeatures, axis=1)

        # Find the index at which the matched vertex would lie, after sorting
        sortedDists = np.argsort(dists)
        sortedInd = np.where(sortedDists == indB)[0][0]
        relativeInd = float(sortedInd) / nRows

        result.append(relativeInd)

    return result


# Performs all appropriate pairwise calls to evaluateAccuracy and returns
# the merged result
def evaluateAllAccuracies(pairList, featureDirectory, accuracyDirectory, correspondenceDirectory = ""):


    skipCount = 0
    for meshA,meshB in pairList:
           
        meshAFilename = os.path.join(featureDirectory, meshA + ".features")
        meshBFilename = os.path.join(featureDirectory, meshB + ".features")

        # The path to put the result at
        resultFilename = os.path.join(accuracyDirectory, meshA + "-" + meshB + ".accuracy")
        
        # Skip computation if the file exists 
        if(os.path.exists(resultFilename)):
            if(args.force):
                os.remove(resultFilename)
            else:
                # print("Skipping accuracy computation for " + meshA + "---" + meshB + ", file already present")
                skipCount += 1
                continue

        print("Evaluating pair: " + meshA + ' -- ' + meshB)


        # Assess the result
        # NOTE THAT WE ARE SUBSAMPLING BY A LOT RIGHT NOW
        subsampleRate = 0.001

        if(correspondenceDirectory == ""):
            pairwiseResult = evaluateAccuracy(parseFeatureFile(meshAFilename), parseFeatureFile(meshBFilename), subsampleRate)
        else:
            correspondenceFilename = os.path.join(correspondenceDirectory, meshA + "-" + meshB + ".groundtruth")
            pairwiseResult = evaluateAccuracyWithFile(parseFeatureFile(meshAFilename), parseFeatureFile(meshBFilename), correspondenceFilename, subsampleRate)

        # Compute precision-recall samples 
        pairwisePR = computePrecisionRecall(pairwiseResult)
   
        # Write the samples to file
        np.savetxt(resultFilename, pairwisePR[1], fmt="%.8f")
                
    if(skipCount > 0):
        print("Skipped evaluating accuracies for " + str(skipCount) + " methods that were already present")

def parseAccuracies(pairList, accuracyDirectory):
  
    print("Reading all accuracies from  " + accuracyDirectory)

    allAccuracies= []
    
    for meshA,meshB in pairList:
   
        # print("Parsing accuracy for pair: " + meshA + ' -- ' + meshB)

        resultFilename = os.path.join(accuracyDirectory, meshA + "-" + meshB + ".accuracy")
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

def plotMethodAccuracy(datasetName, methodName, accuracyLists, plotDir):

    print("Plotting accuracy for method " + methodName)

    yData = np.array(accuracyLists)
    xData = np.linspace(0, 1, num=nPRIncrements)

    sns.set(font_scale=2) 
    plt.title(datasetName + " -- " + methodName + " accuracy")
    plt.xlabel("Threshold")
    plt.ylabel("Fraction under threshold")

    sns.tsplot(yData, xData, err_style="unit_traces", err_palette=sns.color_palette(["#ccd8ff"]))
    # sns.tsplot(yData, xData, ci = [68,95,99.99])
    # sns.tsplot(yData, xData, err_style="boot_traces", n_boot=500)

    outname = plotDir + datasetName + "-" + methodName + "-accuracy.pdf"
    plt.savefig(outname)
    plt.clf()


def plotMerged(averagedData, plotDir):

    print("Plotting merged")

    yData = np.array(averagedData)
    xData = np.linspace(0, 1, num=nPRIncrements)

    sns.set(font_scale=2)
    plt.title("averaged accuracy")
    plt.xlabel("Threshold")
    plt.ylabel("Fraction under threshold")
    
    i = 0
    colors = ["red", "blue", "green", "white", "black"]
    for method in methods:
        sns.tsplot(yData[i], xData, color=colors[i], err_style="unit_traces",
                   err_palette=sns.color_palette(["#ccd8ff"]))
        i = i+1
    
    outname = plotDir + "averaged-accuracy.pdf"
    plt.savefig(outname)
    plt.clf()

### Helpers
def prettyPrintTime(elapsed):

    hours, rem = divmod(elapsed, 3600)
    minutes, seconds = divmod(rem, 60)

    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds)



def parsePairFile(filename):

    print("Parsing pair file " + filename)

    pairs = []

    with open(filename) as pairFile:
        for line in pairFile:

            if line[0] == "#": continue

            items = line.strip().split(",")
            pairs.append((items[0], items[1]))

    print("Found " + str(len(pairs)) + " pairs")

    return pairs 


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
    datasets = ["TOSCA","FAUST"]
    # datasets = ["TOSCA"]

    correspondenceDir = {"TOSCA" : "",
                         "FAUST" : "ground_truth"}

    # Build method lists
    # methods = ["strawman", "curvature", "laplacian-eigenvecs", "hks", "wks"]
    methods = ["strawman", "curvature"]
    methodFunctions = {"strawman" : eval_strawman,
                       "curvature" : eval_curvature,
                       "laplacian-eigenvecs" : eval_laplacian}#,
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
                timingsFilename = os.path.join(outDir, name + ".timing")
                if(os.path.exists(timingsFilename)):
                    if(args.force):
                        os.remove(timingsFilename)
                    else:
                        # print("Skipping method: " + method + " on dataset " + dataset + ". Timings file already present")
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

                timingsFile.write(str(elapsedTime) + "\n")
                timingsFile.close()


    #### Evaluate accuracy ####

    # TODO This code assumes that the ground truth correspondence has
    # the implicit form of equally-numbered verticies matching up. This is
    # just one possible way for a correspondence to be represented, and should
    # be generalized.

    for dataset in datasets:

        # Handle directories for this dataset and method
        datasetPath = os.path.join(dataRoot, dataset)

        # Parse the file that tells us how the dataset is partitioned
        # in to pairs
        pairFilename = os.path.join(datasetPath, "pairs.txt")
        pairList = parsePairFile(pairFilename)

        # Evaluate each method
        for method in methods:

            featuresDir = os.path.join(featuresRoot, dataset, method)
            evaluteDir = os.path.join(evaluateRoot, dataset, method)
            
            if(correspondenceDir[dataset] == ""):
                correspondenceDirectory = ""
            else:
                correspondenceDirectory = os.path.join(datasetPath, correspondenceDir[dataset])


            # Create output directory as needed
            if(not os.path.exists(evaluteDir)):
                os.makedirs(evaluteDir)


            evaluateAllAccuracies(pairList, featuresDir, evaluteDir, correspondenceDirectory)




    #### Generate plots ####
   
    # Make an accuracy plot for each dataset/method
    for dataset in datasets:

        # Handle directories for this dataset and method
        datasetPath = os.path.join(dataRoot, dataset)

        # Parse the file that tells us how the dataset is partitioned
        # in to corresponding pairs
        pairFilename = os.path.join(datasetPath, "pairs.txt")
        pairList = parsePairFile(pairFilename)

        methodAccuracies = []

        for method in methods:

            evaluteDir = os.path.join(evaluateRoot, dataset, method)
   
            # Read in the evaluation results     
            evalResults = parseAccuracies(pairList, evaluteDir)

            # Make a plot of the evaluation results
            plotMethodAccuracy(dataset, method, evalResults, plotRoot)
            A = np.array(evalResults)
            A = np.sum(A, axis=0) / float(len(evalResults))
            methodAccuracies.append(A.tolist())

        # Accumulate to make a merged plot
        plotMerged(methodAccuracies, plotRoot)

        # TODO ROHAN: Performance vs accuracy - area under curve


if __name__ == "__main__":
    main()
