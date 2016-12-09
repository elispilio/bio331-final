import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import csv
import numpy as np
import pylab
import networkx as nx

"""
weights reader and interactome reader
weights reader takes the table file that has the most expressed genes and turns it into a list of the gene names, and the interactome reader takes in and formats the interactome data.
"""
def main():
	edges = reader('interactome.csv')
	nxgraph = nx.Graph()
	nxgraph.add_edges_from(edges)
	nodes = nxgraph.nodes()
	weighted = w_reader('table-4.csv')
	del weighted[0]
	weights = []
	for a in weighted:
		weights.append(a[0])
	graph = [nodes, edges]
	#a = cliqueFinder(graph,weights)
	#list(nx.find_cliques(nxgraph))
	#print(a)
	

def w_reader(filename):
	with open(filename, newline='') as csvfile:
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
	return realdata

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

"""
clique_finder by Eli Spiliotopoulos
the purpose of this set of functions is to read in interactome data and find cliques based on expression, which is expressed by the weights attributed to the nodes. These weights are found via expression data. In the absence of the full RNA-seq dataset to find expression rates, the table of most expressed protiens will be substituted. 
"""
def cliqueFinder(G,weights):
	"""
	finds cliques based on read in weights file
	weight of one attributed to nodes for each hit, in this case of only a few high weight nodes, all others will be 0. Weights is a list of the most expressed nodes. They'll be added to the top of the degree list so cliques will be found using those nodes first, instead of the highest degree node as the normal clique finding algorithm does.
	"""
	nodes = G[0]
	edges = G[1]
	unseen = []
	for n in weights:
		for y in nodes:
			if n == y:
				unseen.append(n)
	deg = {}
	for n in nodes:
		i = 0
		for e in edges:
			if e[0] == n:
				i += 1
			if e[1] == n:
				i += 1
		deg[n] = i

	while True:
		m = 0
		maxi = 0
		for n in nodes:
			if deg[n] > 0 and n not in unseen:
				if deg[n] > m:
					m = deg[n]
					maxi = n
		
		if maxi != 0:
			unseen.append(maxi)
		else:
			break
	neighbs = {}
	for n in nodes:
		neigh = []
		for e in edges:
			if e[0] == n:
				neigh.append(e[1])
			if e[1] == n:
				neigh.append(e[0])
		neighbs[n]=neigh
	setList = []
	while unseen:
		queue = [unseen[0]]
		clusterSet = set()
		while queue:
			for y in queue:
				if y in unseen:
					unseen.remove(y)
			viables = set()
			clusterSet.add(queue[0])
			for n in neighbs[queue[0]]:
				if n in unseen:
					viables.add(n)
			maxim = 0
			m = 0
			neighbViable = set()
			for a in viables:
				if a in unseen:
					if deg[a] > maxim:
						maxim = deg[a]
						m = a

			if m != 0:
				for n in neighbs[m]:
					if n in unseen:
						neighbViable.add(n)
			itera = viables.intersection(neighbViable)
			iterab = []
			for z in itera:
				iterab.append(z)
			for w in iterab:
				queue.append(w)
			del queue[0]
		setList.append(clusterSet)
	return clusterSet

if __name__ == "__main__":
	main()
