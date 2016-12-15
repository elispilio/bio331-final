import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pylab
from itertools import combinations
import networkx as nx
import csv
"""
a simple function to make  images of cliques. also creates associated file of gene names and other ID's
"""
def main():
	d = cliquereader('reportcliques.txt')
	newlist = []
	for cliqueset in d:
		if len(cliqueset) > 2 and len(cliqueset) < 300:
			newlist.append(cliqueset)
	i = 1
	while newlist:
		newgraph = nx.Graph()
		cliques = []
		b = newlist[0]
		for y in combinations(b,2):
			cliques.append(y)
		newgraph.add_edges_from(cliques)
		labels = {}
		for node in newgraph.nodes():
			labels[node]=node
		pos = nx.spring_layout(newgraph)
		for node in newgraph.nodes():
			nx.draw_networkx_nodes(newgraph,pos,nodelist=list(newgraph.nodes()),node_color = '#FF0000', node_size = 100)
			nx.draw_networkx_edges(newgraph,pos,edgelist=list(newgraph.edges()),width=1.0,alpha=0.5,edge_color='#000000')
			nx.draw_networkx_labels(newgraph,pos,labels,font_size=16)
			plt.savefig(str(i)+'.png')
		i+=1
		del newlist[0]

	
def cliquereader(filename):
	with open(filename) as f:
		csvreader = csv.reader(f)
		nodelist = []
		for line in csvreader:
			nodelist.append(line)
	return nodelist

def reader(filenames):

	with open(filenames, newline='') as csvfile:
		spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
		mydata=list(spamreader)
	realdata = []
	for sublist in mydata:
		adder = []
		for item in sublist:
			if item[:7] == 'Unigene':
				adder.append(item)
		if adder != []:
			realdata.append(adder)
	edgelist = []
	for edge in realdata:
		edgelist.append((edge[0], edge[1]))
	
	return edgelist


if __name__ == "__main__":
	main()
