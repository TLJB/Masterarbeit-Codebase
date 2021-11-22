import sys
import re
import os
import getopt
from typing import TextIO
import numpy as np
import math


def read_csv(fobj: TextIO, force: list, disp: list) -> [list, list]:

    nodal_disp = []
    nodal_force = []

    line = fobj.readline()
    line_split = line.split(",")
    for i, name in enumerate(line_split):
        if "DSPL" in name and ":1" in name:
            index_disp = i
        elif "RHS" in name and ":1" in name:
            index_force = i
    line = fobj.readline()
    while line:
        line_split = line.split(",")
        nodal_disp.append(float(line_split[index_disp]))
        nodal_force.append(float(line_split[index_force]))
        line = fobj.readline()

    disp.append(sum(nodal_disp) / len(nodal_disp))
    force.append(sum(nodal_force))

    return [force, disp]


def write_csv(fobj: TextIO, force: list, disp: list):
    fobj.write("Disp Force\n")
    for i in range(len(force)):
        fobj.write("{} {}\n".format(disp[i], abs(force[i])))


def main(argv):

    n_timesteps = 1
    try:
        opts, args = getopt.getopt(argv, "hi:o:n:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print(os.path.basename(__file__) + " -i <inputfile> -o <outputfile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print(
                os.path.basename(__file__) + " -i <inputfile> -o <outputfile>"
            )
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg.strip()
        elif opt in ("-o", "--ofile"):
            outputfile = arg.strip()
        elif opt in ("-n"):
            n_timesteps = int(arg.strip())
    if not os.path.exists(inputfile):
        tb = sys.exc_info()[2]
        raise OSError(
            "inputfile {} not found".format(inputfile)
        ).with_traceback(tb)

    timestep = 0
    force = []
    disp = []
    while timestep < n_timesteps:
        pass
        cfile = inputfile.split("_")[0]
        cfile = cfile + "_" + str(timestep) + ".csv"
        fin = open(cfile)
        force, disp = read_csv(fin, force, disp)
        fin.close()

        timestep += 1

    if not (len(force) == len(disp)):
        tb = sys.exc_info()[2]
        raise ValueError(
            "Length of Force Vector not equal Disp Vector"
        ).with_traceback(tb)

    fout = open(outputfile, "w")
    write_csv(fout, force, disp)
    fout.close()


if __name__ == "__main__":
    main(sys.argv[1:])
