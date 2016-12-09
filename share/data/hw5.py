import operator
#the interactome i was gonna test this on was going to be the GST pull down, which is an experimental technique which finds interactions between known protiens and unknown protiens. It works by becoming expressed if there is an interaction with another protien, which can then be measured via gas chromatography or nuclear magnetic resonance. As I wasnt able to hunt down the bugs in my program I havent properly found a clique example.
def fileReader(filename):
	nodes = []
	edges = []
	with open(filename) as fin:
		for line in fin:
			row = line.strip().split()
			edges.append(row)
	for e in edges:
		if e[0] not in nodes:
			nodes.append(e[0])
		if e[1] not in nodes:
			nodes.append(e[0])

	return nodes,edges

nodes = [1,2,3,4,5,6,7,8,9,10,11,12]
edges = [[1,2],[1,3],[3,1],[3,2],[3,5],[2,5],[5,7],[5,8],[3,4],[4,6],[6,7],[7,8],[7,10],[7,9],[8,10],[8,9],[9,11],[11,12],[12,10],[9,10]]
def greedyPartition(nodes,edges):
	#first make a list of nodes ordered by degree
	#Make a list of unseen nodes
	#loop with a iterator
	#make a cluster list
	#take first node from degree ordered list
	#go through edges and find all neighbors of picked node that are unseen
	#forloop which makes a neighbor list ordered by degree
	#pick first node on that list that is unseen and end inner loop
	#once inner loop is done add that list to a dictionary
	deg = {}
	unseen = []
	for n in nodes:
		i = 0
		for e in edges:
			if e[0] == n:
				i += 1
			if e[1] == n:
				i += 1
		deg[n]=i
	print(deg)
	while len(unseen)< len(nodes):
		#make a list of unseen nodes ordered by degree
		m = 0
		maxi = 0
		for n in nodes:
			if n not in unseen:
				if deg[n] > maxi:
					maxi = deg[n]
					m = n
		unseen.append(m)
	neighbs = {}
	for n in nodes:
		#make a dictionary of neighbors
		neigh = []
		for e in edges:
			if e[0] == n:
				neigh.append(e[1])
			if e[1] == n:
				neigh.append(e[0])
		neighbs[n]=neigh
	"""
	vague pseudocode for the following function
	while no nodes not in sets
		queue= maxdegreenode
		while q:
			viables is cleared
			q[0] is added to cluster set
			take q[0] and add its unseen neighbros to the viable set
			take the max degree neighbor and disjoint neighbViables and viable
			remainder is added, orderd by degree, to queue
			queue[0] is deleted
		clusterset is added LIST
	"""
	setList = []
	while unseen:
		#first while loop, checks if there are nodes that are unseen
		queue = [unseen[0]]
		clusterSet = set()
		while queue:
			#second while loop, takes the maximum degree unseen node, adds its unseen neighbors to a set and compares which other nodes are connected to the nodes in question
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
			#bug is somewhere in here or the next for loop
			for a in viables:
				if a in unseen:
					if deg[a] > maxim:
						maxim = deg[a]
						m = a
			print(m)
			for n in neighbs[m]:
				if n in unseen:
					neighbViable.add(n)
			itera = viables.intersection(neighbViable)
			iterab = []
			for a in itera:
				iterab.append(a)
			for z in iterab:
				maxim = 0
				m = 0
				if deg[z] > maxim:
						maxim = deg[z]
						m = z
				queue.append(m)
				iterab.remove(m)
			del queue[0]
		setList.append(clusterSet)
	print(setList)
