# General purpose harness to run and evaluate algorithms for shape descriptors


import argparse


def main():

    # ==== Argument parsing ====
    parser = argparse.ArgumentParser(description="General harness to run and evaluate shape descriptor algorithms")

    parser.add_argument("root", help="Root path containing data and results")

    parser.add_argument("-f", "--force", help="Overwrite existing results", action="store_true")

    args = parser.parse_args()

    # File location things
    root = args.root;
    dataroot = root + "/data/"

    # Build dataset lists
    datasets = ["TOSCA","FAUST"]
    datasetPaths = {"TOSCA" : dataroot + "tosca/",
                    "FAUST" : dataroot + "faust/"}

    # Build method lists


    # Three stages:
    # 1) Run and benchmark algorithms to generate features on the datasets
    # 2) Evaluate the goodness of each algorithm using the ground truth data
    # 3) Generate all kinds of plots


    #### Run algorithms ####




    #### Evaluate accuracy ####




    #### Generate plots ####




if __name__ == "__main__":
    main()
