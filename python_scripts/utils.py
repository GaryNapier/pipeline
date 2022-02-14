# My general bioinformatics functions

import csv
from csv import DictReader
import re

def flat_list(mylist):
    return [item for sublist in mylist for item in sublist]

def invert_dict_listed(d):
    # Dictionary:
    # x = {'a':['apples'], 'b':['bananas'], 'c':['pears']}
    # invert_dict_listed(x)
    # {'apples': ['a'], 'bananas': ['b'], 'pears': ['c']}

    # x = {'a':['apples'], 'b':['bananas'], 'c':['pears', 'grapes']}
    # invert_dict_listed(x)
    # invert_dict_listed(x)

    inverse = dict()
    for key in d:
        # Go through the list that is saved in the dict:
        for item in d[key]:
            # Check if in the inverted dict the key exists
            if item not in inverse:
                # If not create a new list
                inverse[item] = [key]
            else:
                inverse[item].append(key)
    return inverse

def invert_dict_not_listed(d):
    # Dictionary:
    # x = {'a':'apples', 'b':'bananas', 'c':'pears'}
    # invert_dict_not_listed(x)
    # {'apples': ['a'], 'bananas': ['b'], 'pears': ['c']}

    # Convert to list in case it isn't
    for key in d:
        d[key] = [d[key]]

    inverse = dict()
    for key in d:
        # Go through the list that is saved in the dict:
        for item in d[key]:
            # Check if in the inverted dict the key exists
            if item not in inverse:
                # If not create a new list
                inverse[item] = [key]
            else:
                inverse[item].append(key)
    return inverse

def csv_to_dict(file):

    # If csv has unique rows for each entry of interest (e.g. sample), use this function

    # e.g. 

    # sample, lineage, DR_status
    # S1, lin1, XDR
    # S2, lin3, MDR
    # S3, lin4, sensitive
    # ->
    # {'S1': {'sample':'S1', 'lineage':'lin1', 'DR_status':'XDR'}, 
    # 'S2': {'sample':'S2', 'lineage':'lin3', 'DR_status':'MDR'},
    # 'S3': {'sample':'S3', 'lineage':'lin4', 'DR_status':'sensitive'}}

    reader = csv.DictReader(file)
    key = reader.fieldnames[0]
    out_dict = {}
    for row in reader:
        # Make the id the key, but also put the id in the key-values by including everything
        out_dict[row[key]] = row
    return out_dict

def csv_to_dict_multi(file):
    # If the csv file contains multiple rows for the key of interest ('Gene' in the example below), use this function
    
    # Convert from e.g.:

    # Gene,Mutation,Drug,Confers,Interaction,Literature
    # gyrB,p.Glu540Asp,moxifloxacin,resistance,,10.1128/AAC.00825-17;10.1128/JCM.06860-11
    # pncA,p.Thr100Ile,pyrazinamide,resistance,,10.1128/JCM.01214-17
    # pncA,p.Thr160Ala,pyrazinamide,resistance,,10.1128/JCM.01214-17
    # gid,p.Gly73Ala,streptomycin,resistance,,10.1128/AAC.02175-18
    # gid,p.Leu79Ser,streptomycin,resistance,,10.1128/AAC.02175-18
    # gid,p.Leu108Arg,streptomycin,resistance,,10.1128/AAC.02175-18

    # to:

    # 'rpoC': [{'Gene': 'rpoC',
    #    'Mutation': 'p.Asp485Asn',
    #    'Drug': 'rifampicin',
    #    'Confers': 'resistance',
    #    'Interaction': '',
    #    'Literature': ''},
    #   {'Gene': 'rpoC',
    #    'Mutation': 'p.Asp735Asn',
    #    'Drug': 'rifampicin',
    #    'Confers': 'resistance',
    #    'Interaction': '',
    #    'Literature': ''},...
    # 'rpsL': [{'Gene': 'rpsL',
    #    'Mutation': 'p.Arg86Pro',
    #    'Drug': 'streptomycin',
    #    'Confers': 'resistance',
    #    'Interaction': '',
    #    'Literature': ''},
    #   {'Gene': 'rpsL',
    #    'Mutation': 'p.Arg86Trp',
    #    'Drug': 'streptomycin',
    #    'Confers': 'resistance',
    #    'Interaction': '',
    #    'Literature': ''},... etc

    reader = csv.DictReader(file)
    key = reader.fieldnames[0]
    out_dict = {}
    for row in reader:
        if row[key] in out_dict.keys():
            out_dict[row[key]].append(row)
        else:
            out_dict[row[key]] = []
            out_dict[row[key]].append(row)
    return out_dict

def is_dict(in_dict):
    if not isinstance(in_dict, dict):
        raise Exception('is_dict() function: input is not a dictionary') 

def reformat_mutations(x):
    aa_short2long = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln',
    'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
    'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp',
    'Y': 'Tyr', 'V': 'Val', '*': '*', '-': '-'
    }
    re_obj = re.search("([0-9]+)([A-Z\*])>([0-9]+)([A-Z\*])",x)
    if re_obj:
        codon_num = int(re_obj.group(1))
        ref = aa_short2long[re_obj.group(2)]
        alt = aa_short2long[re_obj.group(4)]
        return "p.%s%s%s" % (ref,codon_num,alt)
    else:
        return None

def resolve_lineages(data):
    # Examples:

    # x = {'lineage2': 5, 'lineage2.2.1':5, 'lineage2.1':5, 'lineage4': 5, 'lineage4.1':5, 'lineage4.1.1':5, 'lineage4.2.1.1':5, 'lineage4.3.1': 5}
    # resolve_lineages(x)
    # {'lineage2': 10,
    # 'lineage2.2.1': 5,
    # 'lineage4': 15,
    # 'lineage4.2.1.1': 5,
    # 'lineage4.3.1': 5}
    
    # x = {'lineage2': 5, 'lineage2.2.1':5, 'lineage2.1':5, 'lineage4': 5, 'lineage4.1':5, 'lineage4.1.1':5, 'lineage4.2.1':5, 'lineage4.3.1':5, 'lineage4.3':5}
    # resolve_lineages(x)
    # {'lineage2': 10, 'lineage2.2.1': 5, 'lineage4': 25, 'lineage4.2.1': 5}

    data = dict(data)
    while True:
        improved = False
        for child in sorted(data,reverse= True,key=lambda x : len(x)):
            parent = ".".join(child.split(".")[:-1])
            if any([parent == x for x in data]):
                data[parent] = data[parent] + data[child] 
                del data[child]
                improved = True
                break
        if not improved:
            break
    return data


def get_embedded_keys(in_dict):
    # Pull keys from embedded dictionary
    
    # Example
    # x = {'samp_1': {'id': 'samp_1', 'lin': 'lin1', 'type': 'XDR'}, 
    # 'samp_2': {'id': 'samp_2', 'lin': 'lin3', 'type': 'MDR'}}
    # x[list(x)[0]].keys()
    # dict_keys(['id', 'lin', 'type'])

    return in_dict[list(in_dict)[0]].keys()