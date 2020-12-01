import argparse
import glob
import math
from os import listdir
from os.path import isfile, join

import matplotlib.pyplot as plt
import scipy.stats as stats
from Levenshtein import distance
from pyLZJD import sim, digest
from pyjarowinkler import distance as jaro_distance

GHIDRA_PATH = ""
R2DEC_PATH = ""
SRC = ""
EXAMPLE_FILE = "alloc.c"
EXAMPLE_FILE_DIR = "C:\\Users\\student\\PycharmProjects\\RadareQUI\\Specimens\\LZJD value incorrect"
BASELINE_DIR = "C:\\Users\\student\\PycharmProjects\\RadareQUI\\Experiment files"
BASELINE_SRC = join(BASELINE_DIR, "BleachedSRC")
BASELINE_OUT = join(BASELINE_DIR, "BleachedOUT")
GHIDRA_NAME = "{}_ghidra_output.txt"
R2DEC_NAME = "{}_r2dec_output.txt"

TOTAL_FILES = 51
# Structures to hold the file values in the form:
DIGESTS = {}
SCORES = {}
WINNERS = {}
KWSEQS_DIGESTS = {}
KWSEQS_SCORES = {}
LEV_SCORES = {}
JARO_SCORES = {}
CKEYWORDS = ['auto', 'break', 'case', 'char', 'const', 'continue', 'default', 'do', 'bool',
             'double', 'else', 'enum', 'extern', 'float', 'for', 'goto', 'if', 'int', 'int64_t',
             'long', 'register', 'return', 'short', 'signed', 'sizeof', 'static', 'struct',
             'switch', 'union', 'unsigned', 'void', 'volatile', 'while']
# change int_64t to just int
# run levenshtein on kw seq

PRINT_SEQ = ['bf_test.c'] #, 'tabulate.c', 'stdfn.c', 'voxelgrid.c', 'vplot.c', 'jitter.c']
        #['plot2d.c', 'save.c']
    # ['external.c']
    # ['internal.c', ]


#note: levenshtein - divide the leven score by the len of the sequence from the original (normalize)

#store kword sequences filter out small ones (under)

#len of kw seq for ghidra in acending other

# do the scatter plot of the sequences levenshtein score
# remove files that are shorter than 10k in size

#chapt 3 is where the data came from and explain the thought procces

#chapt 4 is discusing the results

#look for a levenshtein that factors in the lengths of the strings being compared


#ADD A simple short example (extract out the smallest file)
#Add a discussion on an example where R2 does a better Job than Ghidra and vice versa (show keyword sequences)



def get_keyword_sequence(file):
    #print(" Getting keyword squence for {}".format(file))
    sequence = ""
    with open(file, 'r', encoding='utf-8') as src:
        for line in src:
            for word in line.split():
                if word in CKEYWORDS:
                    sequence += " "
                    if word == "int64_t":
                        sequence += "int"
                        continue
                    sequence += word
    return sequence


def get_args():
    parser = argparse.ArgumentParser(description="Analytics using LZJD on files")
    parser.add_argument('-b', '--baseline', action='store_true', default='store_false', dest='baseline',
                        help='Run Baseline test only')
    parser.add_argument('-e', '--example', action='store_true', default='store_false', dest='example',
                        help='Obtain an example of the process stages')
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
    if isfile(src) and isfile(decompiled):
        with open(src, 'r', encoding='utf=8') as original:
            with open(decompiled, 'r', encoding='utf=8') as dec_output:
                sf = original.read()
                df = dec_output.read()
                lev_score = distance(sf, df)
    else:
        lev_score = distance(src, decompiled)

    return lev_score


def get_jaro_distance(src, decompiled):
    return jaro_distance.get_jaro_distance(src, decompiled, winkler=False)


def plot_boxplt(data, title):
    fig, ax1 = plt.subplots()
    ax1.set_title(title)
    ax1.boxplot(data)
    plt.show()


def plot_hist(data):
    plt.hist(data, bins=math.ceil(math.sqrt(len(data))))
    plt.show()


def plot_scatter(data, title=""):
    filenames = data.keys()
    x = []
    y = []

    for name in filenames:
        x.append(data[name]['ghidra'])
        y.append(data[name]['r2'])

    plt.scatter(x, y)
    plt.xlabel("Ghidra")
    plt.ylabel("R2")
    for i, name in enumerate(data):
        plt.annotate(name, (x[i], y[i]))
    plt.plot((0, 1), 'r--')
    plt.title(title)
    plt.show()


def run_ttest(data1, data2):
    print("Shapiro Output: {}".format(stats.shapiro(data1)))
    print("Shapiro Output: {}".format(stats.shapiro(data2)))
    print("T-test Output: {}".format(stats.ttest_rel(data1, data2)))


def run_keyword_test():
    print("""
        +++++++++++++++++++++++++++++++
        +++ Performing Keyword Test +++
        +++++++++++++++++++++++++++++++
        """)

    for f in listdir(SRC):
        if isfile(join(SRC, f)):
            # print("Keyword test: Running on {}".format(f))
            # prepare a dictionary with the digests ready to compare
            KWSEQS_DIGESTS[f] = {'src': None, 'r2': None, 'ghidra': None}

            # calculate the digest of the source sequence
            # print("Keyword test: Getting digest...".format(f))
            seq_src = get_keyword_sequence(join(SRC, f))

            KWSEQS_DIGESTS[f]['src'] = digest(seq_src)
            # name adjustment
            f2 = f.replace(".c", ".o")

            # calculate digest of ghidra and r2 keyword sequence outputs
            seq_GDR = get_keyword_sequence(join(GHIDRA_PATH, GHIDRA_NAME.format(f2)))
            KWSEQS_DIGESTS[f]['ghidra'] = digest(seq_GDR)
            seq_R2 = get_keyword_sequence(join(R2DEC_PATH, R2DEC_NAME.format(f2)))
            KWSEQS_DIGESTS[f]['r2'] = digest(seq_R2)

            # print_seq = len(seq_src) > 300 and (len(seq_GDR) < 50 or len(seq_R2) < 50)

            # if f in PRINT_SEQ:
            #     print("Sequence for ", f)
            #     print("SRC:", seq_src)
            #     print("GDR:", seq_GDR)
            #     print("R2:", seq_R2)
            # obtain the similarity from source
            # print("Keyword test: Comparing...".format(f))

            KWSEQS_SCORES[f] = {'ghidra': get_lzjd_sim(KWSEQS_DIGESTS[f]['src'], KWSEQS_DIGESTS[f]['ghidra']),
                                'r2': get_lzjd_sim(KWSEQS_DIGESTS[f]['src'], KWSEQS_DIGESTS[f]['r2']),
                                'x': get_lzjd_sim(KWSEQS_DIGESTS[f]['ghidra'], KWSEQS_DIGESTS[f]['r2'])}
    gidra_doms = 0
    for f in KWSEQS_SCORES:
        print("{0:12}: Scores G:{1:20} R2:{2:20} X:{3:20} D:{4:20}".format(f,
                                                                           KWSEQS_SCORES[f]['ghidra'],
                                                                           KWSEQS_SCORES[f]['r2'],
                                                                           KWSEQS_SCORES[f]['x'],
                                                                           KWSEQS_SCORES[f]['ghidra'] -
                                                                           KWSEQS_SCORES[f]['r2'])
              )

        if KWSEQS_SCORES[f]['ghidra'] > KWSEQS_SCORES[f]['r2']:
            gidra_doms += 1
    print("Ghidra Dominated on {} files".format(gidra_doms))

    plot_scatter(KWSEQS_SCORES, title="LZJD Keyword Sequence Scores")

    # obtian the scores as input data to the plots
    bxplt_data_gd = [score['ghidra'] for score in KWSEQS_SCORES.values()]
    bxplt_data_r2 = [score['r2'] for score in KWSEQS_SCORES.values()]

    # create plots
    # plot_boxplt([bxplt_data_gd, bxplt_data_r2], 'Keyword LZJD: Ghidra vs R2')
    # plot_hist(bxplt_data_gd)
    # plot_hist(bxplt_data_r2)

    print("Performing T-Test on LZJD distance of sequences")
    run_ttest(bxplt_data_gd, bxplt_data_r2)


def run_levenshtein_test():
    src_ghidra_lev_scores = []
    src_r2_lev_scores = []
    ghidra_r2_lev_scores = []

    gidra_doms = 0
    for file in listdir(SRC):
        if isfile(join(SRC, file)):
            f = join(SRC, file)
            print("{} Performing Levenshtein Distance on {}".format(gidra_doms, f))

            # remove the extension
            f2 = file.replace(".c", ".o")

            ghidra_file = join(GHIDRA_PATH, GHIDRA_NAME.format(f2))
            r2dec_file = join(R2DEC_PATH, R2DEC_NAME.format(f2))

            src_ghidra_lev_scores.append(get_lev_distance(f, ghidra_file))
            src_r2_lev_scores.append(get_lev_distance(f, r2dec_file))
            ghidra_r2_lev_scores.append(get_lev_distance(ghidra_file, r2dec_file))

        for i in range(TOTAL_FILES):
            print("For file {} LZJD Ghidra:{} R2:{} DIFF:{} both:{}".format(i, src_ghidra_lev_scores[i],
                                                                            src_r2_lev_scores[i],
                                                                            src_ghidra_lev_scores[i] -
                                                                            src_r2_lev_scores[i],
                                                                            ghidra_r2_lev_scores[i]))

            if src_ghidra_lev_scores[i] < src_r2_lev_scores[i]:
                gidra_doms += 1
        print("Ghidra Dominated on {} files".format(gidra_doms))


def run_main_test():
    """
        Test to compare the contents of the source files against the output files of the decompilers
    :return:
    """

    print("""
            +++++++++++++++++++++++++++++++++++++++++++
            +++ Performing Main LZJD Full File Test +++
            +++++++++++++++++++++++++++++++++++++++++++
                """)
    # iterate over the files in the directory
    for f in listdir(SRC):
        if isfile(join(SRC, f)):
            # prepare a dictionary with the digests ready to compare
            DIGESTS[f] = {'src': None, 'r2': None, 'ghidra': None}

            # calculate digest of src file
            DIGESTS[f]['src'] = digest(join(SRC, f))

            # name adjustment
            f2 = f.replace(".c", ".o")

            # calculate digest of ghidra and r2 outputs
            DIGESTS[f]['ghidra'] = digest(join(GHIDRA_PATH, GHIDRA_NAME.format(f2)))
            DIGESTS[f]['r2'] = digest(join(R2DEC_PATH, R2DEC_NAME.format(f2)))

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
    # This section of code prepares visualizations on the data for easy analysis
    plot_scatter(SCORES, title="LZJD Full File scores")

    # obtian the scores as input data to the plots
    bxplt_data_gd = [score['ghidra'] for score in SCORES.values()]
    bxplt_data_r2 = [score['r2'] for score in SCORES.values()]

    # run pairwise t test
    print("Performing T-Test on LZJD Distance of files")
    run_ttest(bxplt_data_gd, bxplt_data_r2)


def run_jaro_kw_test():
    print("""
        ++++++++++++++++++++++++++++++++++++
        +++ Performing JARO Keyword Test +++
        ++++++++++++++++++++++++++++++++++++
        """)
    for f in listdir(SRC):
        if isfile(join(SRC, f)):
            seq_src = get_keyword_sequence(join(SRC, f))
            f2 = f.replace(".c", ".o")
            seq_GDR = get_keyword_sequence(join(GHIDRA_PATH, GHIDRA_NAME.format(f2)))
            seq_R2 = get_keyword_sequence(join(R2DEC_PATH, R2DEC_NAME.format(f2)))

            JARO_SCORES[f] = {'ghidra': get_jaro_distance(seq_src, seq_GDR),
                         'r2': get_jaro_distance(seq_src, seq_R2),
                         'x': get_jaro_distance(seq_GDR, seq_R2)}

    gidra_doms = 0
    for f in JARO_SCORES:
        print("{0:12}: Scores G:{1:20} R2:{2:20} X:{3:20} D:{4:20}".format(f,
                                                                           JARO_SCORES[f]['ghidra'],
                                                                           JARO_SCORES[f]['r2'],
                                                                           JARO_SCORES[f]['x'],
                                                                           JARO_SCORES[f]['ghidra'] -
                                                                           JARO_SCORES[f]['r2']))
        if JARO_SCORES[f]['ghidra'] > JARO_SCORES[f]['r2']:
            gidra_doms += 1
    print("Ghidra Dominated on {} files".format(gidra_doms))

    # Prepare plots in function instead
    # This section of code prepares visualizations on the data for easy analysis
    plot_scatter(JARO_SCORES, title="Jaro Keyword sequence scores")

    # obtian the scores as input data to the plots
    bxplt_data_gd = [score['ghidra'] for score in JARO_SCORES.values()]
    bxplt_data_r2 = [score['r2'] for score in JARO_SCORES.values()]

    # create plots
    # plot_boxplt([bxplt_data_gd, bxplt_data_r2], 'Jaro Distance: Ghidra vs R2')
    # plot_hist(bxplt_data_gd)
    # plot_hist(bxplt_data_r2)

    # run pairwise t test
    print("Performing T-Test on Jaro Distnace of sequences")
    run_ttest(bxplt_data_gd, bxplt_data_r2)


def run_levenshtein_kw_test():
    print("""
        +++++++++++++++++++++++++++++++++++++++++++
        +++ Performing Levenshtein Keyword Test +++
        +++++++++++++++++++++++++++++++++++++++++++
            """)

    for f in listdir(SRC):
        if isfile(join(SRC, f)):
            seq_src = get_keyword_sequence(join(SRC, f))
            f2 = f.replace(".c", ".o")
            seq_GDR = get_keyword_sequence(join(GHIDRA_PATH, GHIDRA_NAME.format(f2)))
            seq_R2 = get_keyword_sequence(join(R2DEC_PATH, R2DEC_NAME.format(f2)))


            # Get Normalization coefficient that wont throw the scale out
            normalization_coefficient = len(seq_src)

            if normalization_coefficient < len(seq_GDR):
                normalization_coefficient = len(seq_GDR)

            if normalization_coefficient < len(seq_R2):
                normalization_coefficient = len(seq_R2)

            # this is normalized by the lenght of the original sequence
            LEV_SCORES[f] = {'ghidra': 1 - get_lev_distance(seq_src, seq_GDR)/normalization_coefficient,
                         'r2': 1 - get_lev_distance(seq_src, seq_R2)/normalization_coefficient,
                         'x': 1 - get_lev_distance(seq_GDR, seq_R2)/normalization_coefficient}

        # if f in PRINT_SEQ:
        #     print("Sequence for ", f)
        #     print("SRC:", seq_src)
        #     print("GDR:", seq_GDR)
        #     print("R2:", seq_R2)
    gidra_doms = 0
    for f in LEV_SCORES:
        print("{0:12}: Scores G:{1:20} R2:{2:20} X:{3:20} D:{4:20}".format(f,
                                                                           LEV_SCORES[f]['ghidra'],
                                                                           LEV_SCORES[f]['r2'],
                                                                           LEV_SCORES[f]['x'],
                                                                           LEV_SCORES[f]['ghidra'] -
                                                                           LEV_SCORES[f]['r2']))
        if LEV_SCORES[f]['ghidra'] > LEV_SCORES[f]['r2']:
            gidra_doms += 1
    print("Ghidra Dominated on {} files".format(gidra_doms))

    # Prepare plots in function instead
    # This section of code prepares visualizations on the data for easy analysis
    plot_scatter(LEV_SCORES, title="Levenshtein Keyword sequence scores")

    # obtian the scores as input data to the plots
    bxplt_data_gd = [score['ghidra'] for score in LEV_SCORES.values()]
    bxplt_data_r2 = [score['r2'] for score in LEV_SCORES.values()]

    # create plots
    # plot_boxplt([bxplt_data_gd, bxplt_data_r2], 'Levenshtein Distance: Ghidra vs R2')
    # plot_hist(bxplt_data_gd)
    # plot_hist(bxplt_data_r2)

    # run pairwise t test
    print("Performing T-Test on Levenshtein Distnace of sequences")
    run_ttest(bxplt_data_gd, bxplt_data_r2)


def run_example():
    print("""
            ++++++++++++++++++++++++++++++
            +++ Performing Example Run +++
            ++++++++++++++++++++++++++++++
            """)
    f = EXAMPLE_FILE
    EKWSEQS_DIGESTS = {}
    EKWSEQS_SCORES = {}
    if isfile(join(EXAMPLE_FILE_DIR, f)):
        # prepare a dictionary with the digests ready to compare
        EKWSEQS_DIGESTS[f] = {'src': None, 'r2': None, 'ghidra': None}


        print("Performing Keyword Sequence Extraction on {}".format(f))
        seq_src = get_keyword_sequence(join(EXAMPLE_FILE_DIR, f))
        print("Keyword Sequence is:\n {}\n".format(seq_src))

        # calculate the digest of the source sequence
        EKWSEQS_DIGESTS[f]['src'] = digest(seq_src)
        # name adjustment
        f2 = f.replace(".c", ".o")

        # calculate digest of ghidra and r2 keyword sequence outputs
        print("Performing Keyword Sequence Extraction on {}".format(GHIDRA_NAME.format(f2)))
        seq_GDR = get_keyword_sequence(join(EXAMPLE_FILE_DIR, GHIDRA_NAME.format(f2)))
        print("Keyword Sequence is:\n {}\n".format(seq_GDR))

        EKWSEQS_DIGESTS[f]['ghidra'] = digest(seq_GDR)

        print("Performing Keyword Sequence Extraction on {}".format(R2DEC_NAME.format(f2)))
        seq_R2 = get_keyword_sequence(join(EXAMPLE_FILE_DIR, R2DEC_NAME.format(f2)))
        print("Keyword Sequence is:\n {}\n".format(seq_R2))

        EKWSEQS_DIGESTS[f]['r2'] = digest(seq_R2)

        print("Performing Keyword Sequence similarity analysis with LZJD ".format(R2DEC_NAME.format(f2)))
        EKWSEQS_SCORES[f] = {'ghidra': get_lzjd_sim(EKWSEQS_DIGESTS[f]['src'], EKWSEQS_DIGESTS[f]['ghidra']),
                            'r2': get_lzjd_sim(EKWSEQS_DIGESTS[f]['src'], EKWSEQS_DIGESTS[f]['r2']),
                            'x': get_lzjd_sim(EKWSEQS_DIGESTS[f]['ghidra'], EKWSEQS_DIGESTS[f]['r2'])}
        for f in EKWSEQS_SCORES:
            print("{0:12}: Scores G:{1:20} R2:{2:20} X:{3:20} D:{4:20}".format(f,
                                                                               EKWSEQS_SCORES[f]['ghidra'],
                                                                               EKWSEQS_SCORES[f]['r2'],
                                                                               EKWSEQS_SCORES[f]['x'],
                                                                               EKWSEQS_SCORES[f]['ghidra'] -
                                                                               EKWSEQS_SCORES[f]['r2']))


def run_get_winners():
    print("""
            +++++++++++++++++++++++++++++
            +++ Obtaining The Winners +++
            +++++++++++++++++++++++++++++
            """)

    for f in KWSEQS_SCORES:
        WINNERS[f] = {}
        if KWSEQS_SCORES[f]['ghidra'] > KWSEQS_SCORES[f]['r2']:
            WINNERS[f]['LZJD'] = 'Ghidra'
        else:
            WINNERS[f]['LZJD'] = 'R2'
    for f in LEV_SCORES:
        if LEV_SCORES[f]['ghidra'] > LEV_SCORES[f]['r2']:
            WINNERS[f]['Levenshtein'] = 'Ghidra'
        else:
            WINNERS[f]['Levenshtein'] = 'R2'
    for f in JARO_SCORES:
        if JARO_SCORES[f]['ghidra'] > JARO_SCORES[f]['r2']:
            WINNERS[f]['Jaro'] = 'Ghidra'
        else:
            WINNERS[f]['Jaro'] = 'R2'


    print("{0:12} & {1:20} & {2:20} & {3:20} ".format('File', 'LZJD', 'Levenshtein', 'Jaro'))
    for f in WINNERS:

        print("{0:12} & {1:20} & {2:20} & {3:20} ".format(f,
                                                          WINNERS[f]['LZJD'],
                                                          WINNERS[f]['Levenshtein'],
                                                          WINNERS[f]['Jaro']))


def main(args):
    # grabbing parameters form arg parse
    global GHIDRA_PATH
    GHIDRA_PATH = args.ghidra_decompiled_path

    global R2DEC_PATH
    R2DEC_PATH = args.r2dec_decompiled_path

    global SRC
    SRC = args.srcpath

    show_lev = args.show_lev
    baseline = args.baseline
    example = args.example

    if baseline is True:
        #####################
        ### Baseline test ###
        #####################

        # This test runs compares two files. The two files the same (original source code in C). One has been tampered
        # to have function names and variable names changed to something similar to what the Ghidra output would look
        # like. This is to give an idea of how big of a similarity score would we have if the decompiler perfectly
        # decompiled code.

        # get baseline source and bleached digests
        baseline_src_digest = get_lzjd_digest(BASELINE_SRC)
        baseline_bleached_digest = get_lzjd_digest(BASELINE_OUT)

        # compare digests with the similarity function from LZJD and output the result
        print("Baseline test performed: LZJD Score for 'ideal decompilation': {}".format(
            get_lzjd_sim(baseline_src_digest[0], baseline_bleached_digest[0])))
        exit(0)

    if example is True:
        run_example()
        exit(0)



    ######################
    ### Main LZJD Test ###
    ######################

    # This test will take all the files in the folders provided in the arguments and will perform an LZJD test
    # with the source vs the output files
    run_main_test()

    ####################
    ## Prepare Plots ###
    ####################

    # This section of code prepares visualizations on the data for easy analysis

    # obtian the scores as input data to the plots
    # bxplt_data_gd = [score['ghidra'] for score in SCORES.values()]
    # bxplt_data_r2 = [score['r2'] for score in SCORES.values()]

    # create plots
    # plot_boxplt([bxplt_data_gd, bxplt_data_r2], 'LZJD: Ghidra vs R2')
    # plot_hist(bxplt_data_gd)
    # plot_hist(bxplt_data_r2)

    #######################
    ### Pairwise T-test ###
    #######################

    # After preparing the plots its a great idea to run a T-test to see if there is a difference in the means of
    # the data. This essentially tells us if one of the decompilers outperformed the other.
    # print("Performing T-Test on Main Test Scores")
    # run_ttest(bxplt_data_gd, bxplt_data_r2)

    #######################
    ### Keyword Testing ###
    #######################

    # This section of code will extract a sequence of c language keywords from source and output files and
    # run the LZJD algorithms again to compare keword sequences.
    run_keyword_test()

    # LEVENSHTEIN on the full files does not make sense....
    # if show_lev is True:
    #     run_levenshtein_test()

    # This will make the same process of extracting the keyword squences and perform the levenshtein distance
    # analysis on the sequences to see how far are they from each other
    run_levenshtein_kw_test()

    # This is just so to see if the Jaro algorithms coincides with LZJD or the Levenshtein scores
    run_jaro_kw_test()

    run_get_winners()

    print("Done.")


if __name__ == '__main__':
    main(get_args())
