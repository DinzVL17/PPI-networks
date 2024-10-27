'''
Q1: Implementing the majority voting network-based candidate protein prediction algorithm
Author: Dinuri Vishara lokuwalpola
Date : 18/1/2023
'''

import networkx as nx

# list the known proteins
AT_stress_proteins= []
with open ("AT_stress_proteins.txt",'r') as file1:
    for lines in file1:
        lines= lines.strip().split("\t")
        lines = lines
        AT_stress_proteins.append(lines[1].upper())

# construct protein-protein interaction network
with open ("string_interactions.tsv",'r') as file2:
    allProteins= []
    pn = nx.Graph()
    for lines in file2:
        if lines != "\n":
            if "#" not in lines:
                lines = lines.strip().split("\t")
                pn.add_edge(lines[0].upper(), lines[1].upper())
                # make a list of all proteins from the network
                allProteins= list(pn.nodes)


# degree of the ATDREB2A
print("degree of the ATDREB2A:", pn.degree("DREB2A"))

# make a list of unknown proteins
unknown_proteins= []
unknown_proteins = set(allProteins) - set(AT_stress_proteins)

# number of unknown proteins in the network for stress tolerance
print("number of unknown proteins:",len(unknown_proteins))

neighbor_dict={}
for n in pn.nodes:
    if n in unknown_proteins:
        # get all the neighbours of unknown proteins
        neighbour = set(nx.neighbors(pn,n))
        # select only known proteins from neighbour list
        new_list = neighbour.intersection(AT_stress_proteins)
        # create a dictionary with unknown protein and its score
        neighbor_dict[n] = len(new_list)
print(neighbor_dict)

# sort the dictionary in descending order by value
from collections import OrderedDict
sorted_dict = sorted(neighbor_dict.items(), key=lambda kv: kv[1], reverse=True)
print(sorted_dict)

# write the sorted dictionary into a file
f = open("sorted_list.txt",'w')
f.write(str(sorted_dict))
f.close()










