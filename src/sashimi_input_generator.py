import sys
import csv
import argparse




def main():



    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files",  help="Comma separated list of file names", required=True)
    parser.add_argument("-p", "--paths",  help="Comma separated list of file paths", required=True)
    parser.add_argument("-c", "--conditions",  help="Comma separated list of condition names", required=True)
    parser.add_argument("-wd", "--working_directory",  help="full path to the snakemake directory", required=True)



    args = parser.parse_args()

    for file, path, condition  in zip(args.files.split(","), args.paths.split(","), args.conditions.split(",")):
        print( "\t".join([ file, args.working_directory + path, condition]) )

main()
