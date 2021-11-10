"""
Find all targets for each drug in drug protein associatin dataset
for each drug find the first order neighbours for proteins interacting with that drug

Note: even though a single drug targets 2 different proteins, the protein neighbor clusters for these proteins might be different.

"""
# This script converts all the uniprot names to their respective grapg label numbers


#!/usr/bin/python python3

import pandas as pd
import sys, argparse, os
import networkx as nx
from tqdm import tqdm

from ensp_to_uniport_stitch import *

# Step1 : convert ENSP-IDs to UniProt-IDs
final_df = ensp_to_uniport()

# Step2: Find all interacting proteins for a drug

drug_prots = dict() # Format of the dictionary {CID:[Prot1, Prot2, .....], CID2:[Prot1, Port3, ....]}
cids = final_df["CIDs"].unique()
for cid in cids:
	drug_prots[cid] = final_df[final_df["CIDs"] == cid]["ENSP-ID"].values


#Step3: Load the graph
print("reading graph...")
Graph = nx.read_gexf("../subgraph.gexf")

#Step4 : Convert all the uniprots to label
print("converting to lables....")
def convert_uniprot_to_label(Graph, drug_prots):
	nodes = dict(Graph.nodes(data=True))
	nodes_1 = dict()
	for n, v in nodes.items():
		nodes_1[nodes[n]["UniProt-ID"]] = nodes[n]["label"]
	nodes = nodes_1
	del(nodes_1)

	edges = list(Graph.edges(data=True))
	edges1 = dict()
	for n in range(len(edges)):
		edges1[edges[n][0:2]] = round(edges[n][-1]["Score"], 5)

	edges = edges1
	del(edges1)

	drugs = {d:[] for d in list(drug_prots.keys())}

	for drug, targets in drug_prots.items():
		for prot in targets:
			if prot in nodes.keys():
				drugs[drug].append(nodes[prot])	
			elif prot not in nodes.keys():
				drugs[drug].append(prot)
			else:
				break

	del(drug_prots)
	
	return nodes, edges, drugs

print("computing new weights....")
nodes, edges, drugs = convert_uniprot_to_label(Graph,drug_prots)
for source, target, attributes in Graph.edges(data=True):
	attributes["Score"] = round((1.01-attributes["Score"]), 5)

print("running Dijkstra...")
for i, _ in tqdm(enumerate(list(drugs.keys()))):
	for j, _ in enumerate(list(drugs.keys())):
		try:
			if list(drugs.keys())[i] != list(drugs.keys())[j]:
				for target1 in drugs[list(drugs.keys())[i]]:
					for target2 in drugs[list(drugs.keys())[j]]:
						target2_neighbors = Graph.neighbors(target2)
						for n in target2_neighbors:
							length, path = nx.single_source_dijkstra(Graph, target1, n, weight="Score")
							print(length, path)
			#elif list(drugs.keys())[i] == list(drugs.keys())[j]:
			#	break
			else:
				break
		except:
			break

