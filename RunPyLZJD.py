import setuptools
import glob
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pdb
import glob

from os import listdir
from os.path import isfile, join
from pyLZJD import sim, digest
from Levenshtein import distance

GHIDRA_PATH = ""
R2DEC_PATH = ""
SRC = ""

GHIDRA_NAME = "{}_ghidra_output.txt"
R2DEC_NAME = "{}_r2dec_output.txt"

TOTAL_FILES = 51

def get_args():
    parser = argparse.ArgumentParser(description="Analytics using LZJD on files")
    parser.add_argument('srcpath', help='Path to the folder containing the source files')
    parser.add_argument('ghidra_decompiled_path', help='Path to the folder containing the ghidra decompiled output files')
    parser.add_argument('r2dec_decompiled_path', help='Path to the folder containing the r2dec decompiled output files')

    return parser.parse_args()


def get_lzjd_digest(path):
    li = glob.glob(join(path, "*.*"))
    return digest(li)


def get_lzjd_sim(src_hash, decompiled_hash):
    return sim(src_hash, decompiled_hash)


def get_lev_distance(src, decompiled):
    lev_score = -1
    with open(src, 'r') as original:
        with open(decompiled, 'r') as dec_output:
            sf = original.read()
            df = dec_output.read()
            lev_score = distance(sf, df)
    return lev_score


def main(args):
    global GHIDRA_PATH
    GHIDRA_PATH = args.ghidra_decompiled_path

    global R2DEC_PATH
    R2DEC_PATH = args.r2dec_decompiled_path

    global SRC
    SRC = args.srcpath

    #src_files = [exe for exe in listdir(args.srcpath) if (isfile(join(args.srcpath, exe)) and '.c' in exe)]
    #ghidra_files = [exe for exe in listdir(GHIDRA_PATH) if (isfile(join(GHIDRA_PATH, exe)) and 'ghidra' in exe)]
    #r2dec_files = [exe for exe in listdir(R2DEC_PATH) if (isfile(join(R2DEC_PATH, exe)) and 'r2dec' in exe)]


    # get the LZJD Digest values for all files
    src_hashes = get_lzjd_digest(SRC)
    ghidra_hashes = get_lzjd_digest(GHIDRA_PATH)
    r2dec_hashes = get_lzjd_digest(R2DEC_PATH)

    # empty scores
    src_ghidra_lzjd_scores = []
    src_r2_lzjd_scores = []
    ghidra_r2_lzjd_scores = []

    src_ghidra_lev_scores = []
    src_r2_lev_scores = []
    ghidra_r2_lev_scores = []
    gidra_doms =0
    for i in range(TOTAL_FILES):

        src_ghidra_lzjd_scores.append(get_lzjd_sim(src_hashes[i], ghidra_hashes[i]))
        src_r2_lzjd_scores.append(get_lzjd_sim(src_hashes[i], r2dec_hashes[i]))
        ghidra_r2_lzjd_scores.append(get_lzjd_sim(ghidra_hashes[i], r2dec_hashes[i]))
    for i in range(TOTAL_FILES):
        print("For file {} LZJD Ghidra:{} R2:{} DIFF:{} both:{}".format(i, src_ghidra_lzjd_scores[i],
                                                                        src_r2_lzjd_scores[i],
                                                                        src_ghidra_lzjd_scores[i] - src_r2_lzjd_scores[
                                                                            i],
                                                                        ghidra_r2_lzjd_scores[i]))

        if src_ghidra_lzjd_scores[i] - src_r2_lzjd_scores[i] >= 0:
            gidra_doms += 1
    print("Ghidra dominates on {} files".format(gidra_doms))

    for file in listdir(SRC):

        if isfile(join(SRC, file)):
            src_file = join(SRC, file)

            # remove the extension
            file = file.replace(".c", ".o")

            ghidra_file = join(GHIDRA_PATH, GHIDRA_NAME.format(file))
            r2dec_file = join(R2DEC_PATH, R2DEC_NAME.format(file))

            src_ghidra_lev_scores.append(get_lev_distance(src_file, ghidra_file))
            src_r2_lev_scores.append(get_lev_distance(src_file, r2dec_file))
            ghidra_r2_lev_scores.append(get_lev_distance(ghidra_file, r2dec_file))



    print("done")


if __name__ == '__main__':
    main(get_args())
