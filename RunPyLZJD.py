import setuptools
import glob
import argparse
import matplotlib.pyplot as plt
import numpy as np

from os import listdir
from os.path import isfile, join
from pyLZJD import sim, digest
from Levenshtein import distance

def get_args():
    parser = argparse.ArgumentParser(description="Analytics using LZJD on files")
    parser.add_argument('srcpath', help='Path to the folder containing the source files')
    parser.add_argument('ghidra_decompiled_path', help='Path to the folder containing the ghidra decompiled output files')
    parser.add_argument('r2dec_decompiled_path', help='Path to the folder containing the r2dec decompiled output files')

    return parser.parse_args()


def get_lzjd_digest(path):
    return digest(path)


def get_lzjd_sim(src_hash, decompiled_hash):
    return sim(src_hash, decompiled_hash)


def get_lev_distance(src, decompiled):
    lev_score = -1
    with open(src, 'r') as original:
        with open(decompiled, 'r') as dec_output:
            lev_score = distance(original, dec_output)
    return lev_score


def main(args):
    files = [exe for exe in listdir(args.srcpath) if isfile(join(args.path, exe))]

    src_hashes = get_lzjd_digest(args.srcpath)
    ghidra_hashes = get_lzjd_digest(args.ghidra_decompiled_path)
    r2dec_hashes = get_lzjd_digest(args.r2dec_decompiled_path)

    if len(src_hashes) != len(ghidra_hashes) or len(src_hashes) != len(r2dec_hashes) :
        raise Exception

    for i in range(len(src_hashes)):
        print(get_lzjd_sim(src_hashes[i], obj_hashes[i]))

    print("done")


if __name__ == '__main__':
    main(get_args())
