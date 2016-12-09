import graphspace_utils
import json_utils
import random
import math
def main():
	"""
	HW3, completed by Eli Spiliotopoulos. Three graphs shared to hw3 group, shortest distances from lab3, random walk graph, absolute differences. Changing the time makes the difference less variable(more nodes of the same color), and the probability seems to cause some of the more connected nodes to be less deeply colored.
	"""
	egfr_r = open("EGFR1-reachable.txt",'r+')
	egfr = []
	for line in egfr_r:
		egfr.append(line.split())
	###read in the file and list out the interactions
	egfr_r.close()
	###calling functions
	nodes,neighbors,types = formatter(egfr)
	graph = (nodes,neighbors)
	walker = rwr(graph,'EGF',.9,1000000)
	edges = []
	for sub in egfr:
		edges.append([sub[0],sub[1]])
	normalshortest = normalizer(bfs('EGF',(nodes,neighbors)))
	nodeAttrs = getNodeAttributes(nodes,abs_dif(walker,normalshortest))
	edgeAttrs = getEdgeAttributes(neighbors,types)
	data = json_utils.make_json_data(nodes,edges,nodeAttrs,edgeAttrs,'absoluteValueDifference','absoluteValueDifference',['HW3'])
	json_utils.write_json(data,'absValDiff.json')
	graphspace_utils.postGraph('absoluteValueDifference','absValDiff.json','elispilio@reed.edu','shareabledrowssap')
	
	return
def abs_dif(walker,normshort):
	###takes the absolute value of the difference between the two normalized sets of values, walker the random walked with repeats and normshort the normalized distances from the bfs
	connect_names = list(walker.keys())
	new_dict = {}
	for n in connect_names:
		new_dict[n] = abs(walker[n]-normshort[n])
	return new_dict

def formatter(connections):
	"""
	parts of this code are taken from the lab3 solution, and creates a list of nodes, a dictionary of neighboring nodes and a dictionary of interaction types
	"""
	egf = {}
	neighbors = {}
	names = []
	for row in connections:
		if (row[0],row[1]) not in egf and row[0] != row[1]:
			egf[(row[0],row[1])]=[row[2]]
		names.append(row[0])
		names.append(row[1])
		if row[0] not in neighbors and row[0] != row[1]:
			neighbors[row[0]]=[row[1]]
		elif row[0] in neighbors and row[0] != row[1]:
			neighbors[row[0]]+=[row[1]]
			
	Rnames = []
	for name in names:
		if name not in Rnames:
			Rnames.append(name)
	names = Rnames
####MAKE AN ADJACENCY MATRIX DO IT DO IT
#### NVM D00D JUST MAKE A DICTIONARY OF EACH NODES INTERACTIONS LOL
	return names,neighbors,egf
	
def rwr(G,S,Q,T):
	"""
	randomly walks through the graph, restarts when probability of above .9 is reached. also restarts when no outgoing nodes are found. T is times to run 
	through the function, originally seeded to get a better idea of what 
	was going on.
	"""
	#random.seed("morehashablestrings")
	counts ={}
	current_w = []
	returnable = {}
	current = S
	nodes,edges = G
	for i in range(T):
		if current not in counts:
			counts[current] = 1
		elif current in counts:
			counts[current] = counts[current]+1

		if random.random() > Q or current not in edges:
			current = S
		else:
			edge_list = edges[current]
			current = random.choice(edge_list)
	
	for n in counts:
		returnable[n] = math.log(counts[n])
	return normalizer(returnable)

def normalizer(counts):
	"""
	for a dictionary with nodes as keys and values as log(counts), perform
	the function 1-(count-absolutemin/absmax-absmin) on every value
	"""
	keys = list(counts.keys())
	vals = list(counts.values())
	normalized = {}
	for node in keys:
		normalized[node]=1-((counts[node]-min(vals))/(max(vals)-min(vals)))
	
	return normalized
def hexmaker(value):
	###takes a value from 0-1 and outputs a hex value for color, ugly copy from lab3
	r = int(.15*255)
	b = int((1-value)*255)
	g = int(.15*255)
	return '#{:02x}{:02x}{:02x}'.format(r,g,b)

def bfs(source,graph):
	###the third time rewritten breadth first search, runs through and notes each level of distance from the source node
	queue = [source]
	nodes = graph[0]
	edges = graph[1]
	testerlist = []
	curated = {source:0}
	###given a dict of node:[edge,edge]
	counter = 0
	while queue:
		varnode = queue[0]
		counter = curated[queue[0]]
		if varnode not in edges:
			pass
		else:
			edgelist = edges[varnode]
		for e in edgelist:
			if e not in curated:
				curated[e]=counter+1
				queue.append(e)
		del queue[0]
	return curated

def getNodeAttributes(nodes,counts):
	nodeAttrs = {}
	for n in nodes:
		if n not in counts:
			counts[n]=0
		nodeAttrs[n] = {}
		nodeAttrs[n]['id'] = n
		nodeAttrs[n]['content'] = n
		nodeAttrs[n]['text_border_color'] = '#FFFFFF'
		###Idk if i'm inputting it wrong but textbordercolor doesnt work right
		nodeAttrs[n]['height'] = 50
		nodeAttrs[n]['width'] = 50
		nodeAttrs[n]['background_color'] = hexmaker(counts[n]) 
	return nodeAttrs

def getEdgeAttributes(neighbors,egf):
	edges = list(egf.keys())
	attrs = {}
	for e in edges:
		source = e[0]
		target = e[1]
		if source not in attrs:
			attrs[source] = {}
		attrs[source][target] = {}
		attrs[source][target]['popup'] = '<b>'+str(egf[e])+'</b>'
	return attrs

if __name__ == '__main__':

	main()
