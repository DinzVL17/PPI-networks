'''
Title: Implementation of a similarity-based network link prediction algorithm for Drug-Target Interactions (DTI).
Author: Dinuri Vishara
Date: 15/02/2023
'''
import networkx as nx
from itertools import islice

# open the file
with open("DTIsubset.tsv",'r') as file:
    drug_list=[]
    target_list=[]
    pn = nx.Graph()
    # remove the first line of the file
    for lines in islice(file,1,None):
        lines = lines.strip().split("\t")
        if lines != "\n":
            # save drugs into a list
            drug_list.append(lines[0])
            # save the targets into a list
            target_list.append(lines[1])
            # remove duplicates in drugs list
            drugs = set(drug_list)
            # remove duplicates in targets list
            targets = set(target_list)
            # constructing the DTI network
            pn.add_edge(lines[0], lines[1])


diction={}
# create a list of neighbours of drugs
for d in drugs:
    drugs_neigh= list(pn.neighbors(d))

    # create a list of neighbours of targets
    for t in targets:
        # remove already existing interactions
        if t not in drugs_neigh:
            targets_neigh= list(pn.neighbors(t))
            second_neigh_list = []
            # get the neighbours of non-interacting targets and save in a list
            for n in targets_neigh:
                second_neigh = list(pn.neighbors(n))
                second_neigh_list += second_neigh
                # get the common targets of neighbours of drugs and neighbours of drugs of non-interacting targets
                common_neigh = set(second_neigh_list).intersection(set(drugs_neigh))
                # calculate the CN score for each interaction
                score = len(common_neigh)
                # save interactions and CN scores in a dictionary
                diction[t,d]= score
print(diction)

# sort the CN scores for novel drug-target interactions in descending order
sorted_dict = sorted(diction.items(),key=lambda kv:kv[1],reverse=True)
print("CN scores for novel drug-target interactions:\n",sorted_dict)





