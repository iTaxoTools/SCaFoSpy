#### SCaFoSpy V0.1 ####
#### David Leiße ####
#### david.leisse@uni-bielefeld.de ####


import Utils
import sys 
from os import listdir
from os.path import isfile, join, basename

__usage__ = """
            python3 SCaFoS.py
            --dir <PATH_TO_INPUT_DIRECTORY>
            --out <PATH_TO_OUTPUT_DIRECTORY>
            --type <select between 'nucleotide' OR 'protein'>

            """

def compare_sequences_by_character(sequences_list: list, type : str) -> str:
    """
    Comparing multiple sequences by character and creating one consesus sequence.

    :param List sequences_list:  List of sequences to compare of the same length.

    :return consensus: Consensus sequence created from provided sequences
    """
    if type == "protein":
        missing = ["*","?"," "]
    
    if type == "nucleotide":
        missing = ["*","n","N","?"," "]
    consensus = []

    if not sequences_list:
        return print("Sequence list did not contain any sequences!")

    for n in range(len(sequences_list[0])):
        
        na = []
        chars = []
        for seq in sequences_list:
            
            if seq[n] in missing:
                na.append(seq[n])
            
            else:
                chars.append(seq[n])
        
        if len(chars) == 1:
            consensus.append(chars[0])
        
        elif chars:
            if Utils.all_equal(chars):
                consensus.append(chars[0])
            else:
                consensus.append("?")
        
        else:
            consensus.append(Utils.most_frequent(na))

    consensus = "".join(consensus)
    return consensus


def __Main__(args):

    if type(args) == str:
        args = args.split(" ")

    dir = args[args.index("--dir")+1 ]
    outdir = args[args.index("--out") +1 ]
    types = args[args.index("--type")+1 ]

    files = [join(dir, f) for f in listdir(dir) if isfile(join(dir, f))]

    for file in files:
        filename = basename(file)
        aligns = Utils.load_fasta_ali_file(file)

        sorted_aligns = Utils.sort_aligns_to_sample(aligns)

        consensus_seqs = []
        for name,sequences in sorted_aligns.items():
            consensus = compare_sequences_by_character(sequences,types)
            consensus_seqs.append([name,consensus])

        outfile = filename
        Utils.write_scafos_output(outdir,outfile,consensus_seqs,types)

__Main__("--dir /Users/david/Development/SCAFOS/data --out /Users/david/Development/SCAFOS/data --type protein")

"""
if "--dir" in sys.argv and "--out" in sys.argv and "--type" in sys.argv:
    __Main__(sys.argv)

else:
    print(__usage__)
"""