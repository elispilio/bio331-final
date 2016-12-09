import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import csv
import numpy as np
import pylab
import networkx as nx
with open('interactome.csv', newline='') as csvfile:
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
newedges = []

#for x in range(0,100):
#	newedges.append(edgelist[x])
#testest = nx.Graph()
#testest.add_edges_from(newedges)
#print(nx.all_pairs_node_connectivity(testest))
testGraph = nx.Graph()
testGraph.add_edges_from(edgelist)
nodes = testGraph.nodes()
degList = []
for node in nodes:
	degCounter = 0
	for edge in edgelist:
		if edge[0] == node or edge[1] == node:
			degCounter += 1
	degList.append(degCounter)
print(nodes)
#l = degList
#counts = [[x,l.count(x)] for x in set(l)]
#n, bins, patches = pylab.hist(degList,25,normed=1,histtype='stepfilled')
#pylab.setp(patches, 'facecolor', 'g', 'alpha',0.75)
#pylab.savefig('degreedistribution.png', bbox_inches='tight')
#print(nx.average_clustering(testGraph))
