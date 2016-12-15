import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import csv
import numpy as np
import pylab
from itertools import combinations
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
	weighted = w_reader('annotatedweight.txt')
	del weighted[0]
	weights = []
	for a in weighted:
		weights.append(a[0])
	graph = [nodes, edges]
	#p = cliquepics('cliques.txt')
	#newgraph = nx.Graph()
	#for edge in p:
		#if len(edge) == 6:
			#newlist = edge
	#newedges = []
	#for edge in combinations(newlist,2):
		#newedges.append(edge)
	#newgraph.add_edges_from(newedges)
	#labels = {}
	#for node in newlist:
		#labels[node]=node
	ab = cliqueFinder(graph,weights)
	reportcliques = open('reportcliques.txt','w')
	#ab = list(nx.find_cliques(nxgraph))
	for x in ab:
		reportcliques.write(str(x))
		reportcliques.write('\n')
	#newnewnewcliques.close()
	#pos = nx.spring_layout(newgraph)
	#for node in newgraph.nodes():
		#nx.draw_networkx_nodes(newgraph,pos,nodelist=newlist,node_color='#FF0000',node_size = 100, alpha=0.5,with_labels = True)
		#nx.draw_networkx_edges(newgraph,pos,edgelist=newedges,width=1.0,alpha=0.5,edge_color='#000000')
		#nx.draw_networkx_labels(newgraph,pos,labels,font_size=16)
		#plt.savefig('ohgodplswork3.png')
def cliquepics(cliquesfile):
	with open(cliquesfile) as f:
		csvreader = csv.reader(f)
		nodelist = []
		for line in csvreader:
			#l = line.split('""')[1::2]
			nodelist.append(line)
	return nodelist
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
testlist = []
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
	while unseen: #while unseen is not empty
		queue = [unseen[0]] #unseen is ordered by weights and degree first unseen goes to queue
		clusterSet = set() #a new clique is opened
		while queue: #while there are nodes in the queue
			for y in queue:
				if y in unseen:
					unseen.remove(y) #remove the node from unseen
			viables = set()
			clusterSet.add(queue[0]) #the node is the first member of the set
			for n in neighbs[queue[0]]: #all the neighbors that arent seen are viable members of the clique
				if n in unseen:
					viables.add(n)
			maxim = 0
			m = 0
			neighbViable = set()
			testlist.append(viables)
			for a in viables:
				if a in unseen:
					if deg[a] > maxim:
						maxim = deg[a]
						m = a

			if m != 0:
				for n in neighbs[m]:
					if n in unseen:
						neighbViable.add(n) #add the maximum degree neighbor of the viable neighbor
			iterab = list(viables.intersection(neighbViable)) # all the nodes that are connected to both of them
			for w in iterab: #adds to queue by degree
				maxim = 0
				m = 0
				if deg[w] > maxim:
					maxim = deg[w]
					m = w
				queue.append(m) #add the next node to the queue
				iterab.remove(m)#remove the node from iterab
			del queue[0]
		setList.append(clusterSet) #add the clique to the list
	return setList

if __name__ == "__main__":
	main()
