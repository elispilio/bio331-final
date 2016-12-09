"""
Lab6 tsting how to implement the most used functions from the course so far
"""
def fileReader(graphfile):
	graph = open(graphfile,'r+')
	edges = []
	for line in graph:
		edges.append(line.split())
	graph.close
	nodes = []
	for edge in edges:
		if edge[0] not in nodes:
			nodes.append(edge[0])
		elif edge[1] not in nodes:
			nodes.append(edge[1])
	return nodes,edges

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


def frequencyFinder(nodes, edges):
	selfs = []
	ffl = set()
	fbl = set()
	adj = {}
	tester = []
	for node in nodes:
		for edge in edges:
			if edge[0] == node:
				tester.append(edge[1])
		adj[node] = tester
		tester = []
	#self checker
	for node in nodes:
		if node in adj[node]:	
			selfs.append(node)
	#feedback checker
	for node in nodes:
		for neighbor in adj[node]:
			if node not in adj[neighbor]:
				for n2 in adj[neighbor]:
					if neighbor not in adj[n2]:
						if node in adj[n2]:
							if n2 not in adj[node]:
								fbl.add(frozenset([node,neighbor,n2]))
	#lazy feed forward checker
	for node in nodes:
		for neighbor in adj[node]:
			if node not in adj[neighbor]:
				for n2 in adj[neighbor]:
					if neighbor not in adj[node]:
						if n2 in adj[node] and node not in adj[n2]:
							ffl.add(frozenset([node,neighbor,n2]))
	return selfs,fbl,ffl

def hexmaker(value):
	###takes a value from 0-1 and outputs a hex value for color, ugly copy from lab3
	r = int(.15*255)
	b = int((1-value)*255)
	g = int(.15*255)
	return '#{:02x}{:02x}{:02x}'.format(r,g,b)

class graph:
	def __init__(self,nodes,edges):
		self.nodes = nodes
		self.edges = edges

	def get_degree(self,node):
		i = 0
		for edge in edges:
			if node == edge[0]:
				i += 1
		return i

	

	def node_attributes(self,nodes):
		node_attrs = {}
		for node in nodes:
			node_attrs[node] = {}
			node_attrs[node]['degree'] = self.get_degree(node)
			node_attrs[node]['id'] = node
			node_attrs[node]['content'] = node
			node_attrs[node]['border_color'] = 'black'
			node_attrs[node]['height'] = 50
			node_attrs[node]['width'] = 50
		
		return node_attrs

	def edge_attributes(self,edges):
		edge_attrs = {}
		for edge in edges:
			source = edge[0]
			target = edge[1]
			if source not in edge_attrs:
				edge_attrs[source] = {}
				edge_attrs[source][target] = {}
				edge_attrs[source][target]['target_arrow_shape'] = 'triangle'
				edge_attrs[source][target]['target_arrow_fill'] = 'filled'
				edge_attrs[source][target]['target_arrow_color'] = 'black'
		return edge_attrs

nodes,edges = fileReader("testgraph.txt")
aces = graph(nodes,edges)
