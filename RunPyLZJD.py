import setuptools
import glob
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pdb
import glob, math

from os import listdir
from os.path import isfile, join
from pyLZJD import sim, digest
from Levenshtein import distance
import scipy.stats as stats

GHIDRA_PATH = ""
R2DEC_PATH = ""
SRC = ""

BASELINE_DIR = "C:\\Users\\student\\PycharmProjects\\RadareQUI\\Experiment files"
BASELINE_SRC = join(BASELINE_DIR, "BleachedSRC")
BASELINE_OUT = join(BASELINE_DIR, "BleachedOUT")
GHIDRA_NAME = "{}_ghidra_output.txt"
R2DEC_NAME = "{}_r2dec_output.txt"

TOTAL_FILES = 51
# Structures to hold the file values in the form:
#   {'file_name': {'src_ghidra': []  # src-ghidra lzjd values
#                 'src_r2': []  # src-r2 lzjd values
#                 'ghidra_r2': []}  # ghidra-r2 lzjd values }
DIGESTS = {}
SCORES = {}

# Structures to hold the file values in the form:
#   {'file_name': {'src_ghidra': []} # src-ghidra LEV values
#                 {'src_r2': []}  # src-r2 LEV values
#                 {'ghidra_r2': []} # ghidra-r2 LEV values }
FILE_SCORES_LEV = {}


def get_args():
    parser = argparse.ArgumentParser(description="Analytics using LZJD on files")
    parser.add_argument('-b', '--baseline', action='store_true', default='store_false', dest='baseline',
                        help='Run Baseline test only')
    parser.add_argument('-l', '--levenshtein', action='store_true', default='store_false', dest='show_lev',
                        help='Run Levenshtein Distance calculations')

    parser.add_argument('srcpath', help='Path to the folder containing the source files')
    parser.add_argument('ghidra_decompiled_path', help='Path to the folder containing the ghidra decompiled output '
                                                       'files')
    parser.add_argument('r2dec_decompiled_path', help='Path to the folder containing the r2dec decompiled output files')

    return parser.parse_args()


def get_lzjd_digest(path):
    li = glob.glob(join(path, "*.*"))
    return digest(li)


def get_lzjd_sim(src_hash, decompiled_hash):
    return sim(src_hash, decompiled_hash)


def get_lev_distance(src, decompiled):
    lev_score = -1
    with open(src, 'r', encoding='utf=8') as original:
        with open(decompiled, 'r', encoding='utf=8') as dec_output:
            sf = original.read()
            df = dec_output.read()
            lev_score = distance(sf, df)
    return lev_score


def plot_boxplt(data, title):
    fig, ax1 = plt.subplots()
    ax1.set_title(title)
    ax1.boxplot(data)
    plt.show()


def plot_hist(data):
    plt.hist(data, bins=math.ceil(math.sqrt(len(data)) ))
    plt.show()


def run_ttest(data1, data2):
    print("Shapiro Output: {}".format(stats.shapiro(data1)))
    print("Shapiro Output: {}".format(stats.shapiro(data2)))
    print("T-test Output: {}".format(stats.ttest_rel(data1, data2)))


def main(args):
    # grabbing parameters form arg parse
    global GHIDRA_PATH
    GHIDRA_PATH = args.ghidra_decompiled_path

    global R2DEC_PATH
    R2DEC_PATH = args.r2dec_decompiled_path

    global SRC
    SRC = args.srcpath

    show_lev = args.show_lev
    baseline = args.show_lev

    baseline_src_digest = get_lzjd_digest(BASELINE_SRC)
    baseline_bleached_digest = get_lzjd_digest(BASELINE_OUT)
    print("Baseline test performed: LZJD Score for 'ideal decompilation': {}".format(
            get_lzjd_sim(baseline_src_digest[0], baseline_bleached_digest[0])))
    exit(0)

    for f in listdir(SRC):
        if isfile(join(SRC, f)):
            global DIGESTS
            DIGESTS[f] = {'src': None, 'r2': None, 'ghidra': None}

            # calculate digest of src file
            DIGESTS[f]['src'] = digest(join(SRC, f))

            # name adjustment
            f2 = f.replace(".c", ".o")

            # calculate digest of ghidra and r2 outputs
            DIGESTS[f]['ghidra'] = digest(join(GHIDRA_PATH, GHIDRA_NAME.format(f2)))
            DIGESTS[f]['r2'] = digest(join(R2DEC_PATH, R2DEC_NAME.format(f2)))

            global SCORES
            # obtain the similarity from source
            SCORES[f] = {'ghidra': get_lzjd_sim(DIGESTS[f]['src'], DIGESTS[f]['ghidra']),
                         'r2': get_lzjd_sim(DIGESTS[f]['src'], DIGESTS[f]['r2']),
                         'x': get_lzjd_sim(DIGESTS[f]['ghidra'], DIGESTS[f]['r2'])}

    gidra_doms = 0
    for f in SCORES:
        print("{0:12}: Scores G:{1:20} R2:{2:20} X:{3:20} D:{4:20}".format(f,
                                                                SCORES[f]['ghidra'],
                                                                SCORES[f]['r2'],
                                                                SCORES[f]['x'],
                                                                SCORES[f]['ghidra'] - SCORES[f]['r2']))
        if SCORES[f]['ghidra'] > SCORES[f]['r2']:
            gidra_doms += 1
    print("Ghidra Dominated on {} files".format(gidra_doms))

    bxplt_data_gd = [score['ghidra'] for score in SCORES.values()]
    bxplt_data_r1 = [score['r2'] for score in SCORES.values()]

    plot_boxplt([bxplt_data_gd, bxplt_data_r1], 'Ghidra vs R2')
    plot_hist(bxplt_data_gd)
    plot_hist(bxplt_data_r1)

    run_ttest(bxplt_data_gd,bxplt_data_r1)

    if show_lev is True:
        src_ghidra_lev_scores = []
        src_r2_lev_scores = []
        ghidra_r2_lev_scores = []

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

            for i in range(TOTAL_FILES):
                print("For file {} LZJD Ghidra:{} R2:{} DIFF:{} both:{}".format(i, src_ghidra_lev_scores[i],
                                                                                src_r2_lev_scores[i],
                                                                                src_ghidra_lev_scores[i] - src_r2_lev_scores[i],
                                                                                ghidra_r2_lev_scores[i]))

                if src_ghidra_lev_scores[i] < src_r2_lev_scores[i]:
                    gidra_doms += 1
            print("Ghidra Dominated on {} files".format(gidra_doms))

    print("done")


if __name__ == '__main__':
    main(get_args())
