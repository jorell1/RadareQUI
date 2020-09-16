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

TOTAL_FILES = -1
# Structures to hold the file values in the form:
#   {'file_name': {'src_ghidra': []  # src-ghidra lzjd values
#                 'src_r2': []  # src-r2 lzjd values
#                 'ghidra_r2': []}  # ghidra-r2 lzjd values }
DIGESTS = {}

# Structures to hold the file values in the form:
#   {'file_name': {'src_ghidra': []} # src-ghidra LEV values
#                 {'src_r2': []}  # src-r2 LEV values
#                 {'ghidra_r2': []} # ghidra-r2 LEV values }
FILE_SCORES_LEV = {}


def get_args():
    parser = argparse.ArgumentParser(description="Analytics using LZJD on files")
    parser.add_argument('srcpath', help='Path to the folder containing the source files')
    parser.add_argument('ghidra_decompiled_path', help='Path to the folder containing the ghidra decompiled output '
                                                       'files')
    parser.add_argument('r2dec_decompiled_path', help='Path to the folder containing the r2dec decompiled output files')
    parser.add_argument('-l', '--levenshtein', action='store_true', default='store_false', dest='show_lev',
                        help='Run Levenshtein Distance calculations')

    return parser.parse_args()


def get_lzjd_digest(path):
    li = glob.glob(join(path, "*.*"))
    return digest(li)


def get_lzjd_sim(src_hash, decompiled_hash):
    return 1 - sim(src_hash, decompiled_hash)


def get_lev_distance(src, decompiled):
    lev_score = -1
    with open(src, 'r', encoding='utf=8') as original:
        with open(decompiled, 'r', encoding='utf=8') as dec_output:
            sf = original.read()
            df = dec_output.read()
            lev_score = distance(sf, df)
    return lev_score


def prepare_plot():
    figure, ax = plt.subplots()
    ax.plot()


def main(args):
    # grabbing parameters form arg parse
    global GHIDRA_PATH
    GHIDRA_PATH = args.ghidra_decompiled_path

    global R2DEC_PATH
    R2DEC_PATH = args.r2dec_decompiled_path

    global SRC
    SRC = args.srcpath

    show_lev = args.show_lev

    # for f in listdir(SRC):
    #     if isfile(join(SRC, f)):
    #         global DIGESTS
    #         DIGESTS[f] = {'src': None, 'r2': None, 'ghidra': None}
    #
    #         # calculate digest of file
    #         DIGESTS[f]['src'] = digest(join(SRC, f))
    #         f = f.replace(".c", ".o")
    #         DIGESTS[f]['ghidra'] = digest(join(GHIDRA_PATH, GHIDRA_NAME.format(f)))
    #         DIGESTS[f]['r2'] = digest(join(R2DEC_PATH, R2DEC_NAME.format(f)))

    global TOTAL_FILES
    #TOTAL_FILES = len(DIGESTS)

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

    gidra_doms = 0

    for i in range(len):
        src_ghidra_lzjd_scores.append(get_lzjd_sim(src_hashes[i], ghidra_hashes[i]))
        src_r2_lzjd_scores.append(get_lzjd_sim(src_hashes[i], r2dec_hashes[i]))
        ghidra_r2_lzjd_scores.append(get_lzjd_sim(ghidra_hashes[i], r2dec_hashes[i]))

    for i in range(TOTAL_FILES):
        print("For file {} LZJD Ghidra:{} R2:{} DIFF:{} both:{}".format(i, src_ghidra_lzjd_scores[i],
                                                                        src_r2_lzjd_scores[i],
                                                                        src_ghidra_lzjd_scores[i]-src_r2_lzjd_scores[i],
                                                                        ghidra_r2_lzjd_scores[i]))

        if src_ghidra_lzjd_scores[i] - src_r2_lzjd_scores[i] < 0:
            gidra_doms += 1

    print("Ghidra is dominated on {} files".format(gidra_doms))

    gidra_doms = 0

    for file in listdir(SRC):
        if isfile(join(SRC, file)):
            src_file = join(SRC, file)
            print("{} Performing Levenshtein Distance on {}".format(gidra_doms, src_file))

            # remove the extension
            file = file.replace(".c", ".o")

            ghidra_file = join(GHIDRA_PATH, GHIDRA_NAME.format(file))
            r2dec_file = join(R2DEC_PATH, R2DEC_NAME.format(file))

            src_ghidra_lev_scores.append(get_lev_distance(src_file, ghidra_file))
            src_r2_lev_scores.append(get_lev_distance(src_file, r2dec_file))
            ghidra_r2_lev_scores.append(get_lev_distance(ghidra_file, r2dec_file))
            gidra_doms += 1

    gidra_doms = 0
    if show_lev:
        for i in range(TOTAL_FILES):
            print("For file {} LZJD Ghidra:{} R2:{} DIFF:{} both:{}".format(i, src_ghidra_lev_scores[i],
                                                                            src_r2_lev_scores[i],
                                                                            src_ghidra_lev_scores[i] - src_r2_lev_scores[i],
                                                                            ghidra_r2_lev_scores[i]))

            if src_ghidra_lev_scores[i] - src_r2_lev_scores[i] < 0:
                gidra_doms += 1
        print("Ghidra is dominated on {} files".format(gidra_doms))

    prepare_plot()

    plt.plot([i for i in range(TOTAL_FILES)], src_ghidra_lzjd_scores, )
    print("done")


if __name__ == '__main__':
    main(get_args())
