##### Utils.py #####
##### David Leiße #####
##### david.leisse@uni-bielefeld.de #####

import dendropy
from os import mkdir
from genericpath import isdir
import os.path
from itertools import groupby

def load_tree(string: str, mode: str, schema: str) -> dendropy.Tree:
    """
    Loading tree from .tre file and plotting it.
    Working with Newick schema until now.

    :param String string: Relative or absolute path to .tre file or data string
    :param String mode: Choosing between 'path' and 'string'
    :param String schema: I.e. 'newick' or 'nexus' depending on .tre file

    :return tree: Dendropy.Tree object
    """
    if mode == "path":
        tree = dendropy.Tree.get(path = string, schema = schema,rooting ="force-unrooted")
    elif mode == "string":
        tree = dendropy.Tree.get(data = string, schema = schema)
        
    return tree

def load_fasta_ali_file(file: str) -> list:
    """
    Loading fasta file into list. List containing lists with each an ID and sequence

    :param String file: full path to file
    :return ls: List containing lists with each an ID and sequence  
    """

    with open(file, "r") as f:
        line = f.readline()
        characters = ["@"]
        characters2 = [" "]
        sequence = []
        header = ""
        ls = []
        while line:
            
            while line and not line.startswith(">"):
                if not sequence == "\n":
                    sequence.append(line.replace("\n",""))
                line = f.readline()

            if sequence and header:
                ls.append([header, "".join(sequence)])
                
            sequence = []
            header = line[1:].strip()
            
            #### SCAFOS version!
            for a,b in zip(characters,characters2):
                if a in header:
                    header = header.split(a)[0]
                if b in header:
                    header = header.replace(b,"_")
            

            line = f.readline()
    
    return ls

def write_scafos_output(directory: str, file: str, seq_ls: list, type: str):
    """
    Writing output files from modified dictionary

    :param List mod_ali_ls: modified list with contents for decontaminated .ali files
    :param String ali_file: path to original .ali file
    """

    output_directory = directory + "/consensus"
    if not isdir(output_directory):
        mkdir(output_directory)

    file = output_directory + "/" + file 

    with open(file, "w") as out:
        if type == "protein":
            out.write("#\n#\n")
        for id,seq in seq_ls:
            out.write(">" + id + "\n")
            out.write(seq + "\n")

def get_index_in_list(value: any, list: list) -> int:
    """
    Getting index of list of lists by value.

    :param Any value: value that is searched in list of lists
    :param List list: corresponding List of lists

    :return index: Index that corresponds to list that contains value
    """
    index = -1
    for idx,x in enumerate(list):
        if value in x:
            index = idx

    if index == -1:
        return None

    else:
        return index

def get_corresponding_files(tre_files: list, ali_files: list, fasta_files: list) -> list:
    """
    Getting all corresponding files to a .tre file. 
    Corresponding files are name similar to .tre file.
    Usually .fasta / .ali files

    :param List tre_files: List of .tre files
    :param List ali_files: List of .ali files
    :param List fasta_files: List of .fasta files

    :return file_list: List of lists containing .tre and corresponding .fasta and .ali files
    """
    file_list = []
    for tre_file in tre_files:
        temp = []
        name = os.path.splitext(tre_file)[0].replace("-","_")
        temp.append(tre_file)

        for ali_file in ali_files:
            ali = os.path.splitext(ali_file)[0].replace("-","_")
            if ali == name:
                temp.append(ali_file)
                
        for fasta_file in fasta_files:
            fasta = os.path.splitext(fasta_file)[0].replace("-","_")
            if fasta == name:
                temp.append(fasta_file)
        file_list.append(temp)

    return file_list

def sort_aligns_to_sample(alignments: list) -> dict:
    """
    Sorting multiple sequences to one single sample by substring in sequence name.
    For later comparsion of sequences.

    :param List alignments: List of sequence names and their sequences from alignment file (.ali or .fasta)

    :return outdict: Dictionary with sorted sequences. Multiple sequences in a list are sorted to the respective sample name.
    """

    outdict = {}

    for sample in alignments:
        id = sample[0]
        seq = sample[1]
        if not id in outdict:
            outdict.update({id : [seq]})
        else: 
            temp = outdict[id]
            temp.append(seq)
            outdict.update({id : temp})

    return outdict



def all_equal(iterable):
    """
    Checking if in an iterable all items are the same
    """
    g = groupby(iterable)
    return next(g, True) and not next(g, False)

def most_frequent(iterable):
    """
    Checking which item is the most common in an iterable.
    """

    return max(set(iterable), key=iterable.count)