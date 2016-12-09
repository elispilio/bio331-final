"""
Lab 7, implementation of the neighbor joining algorithm
"""
import graphspace_utils
import json_utils
def main():
	D,headers = readFromFile('wikipedia-example.txt')
	tree = neighb_join(D,headers)
	nodeAttrs = getNodeAttributes(tree)
	edgeAttrs = getEdgeAttributes(tree)
	edges = []
	for edge in tree:
		edges.append([edge[0],edge[1]])
	nodes = []
	for node in tree:
		if node[0] not in nodes:
			nodes.append(node[0])
		if node[1] not in nodes:
			nodes.append(node[1])
	data = json_utils.make_json_data(nodes,edges,nodeAttrs,edgeAttrs,'classExampleEli','classExampleEli',['lab7'])
	json_utils.write_json(data,'classExampleLab7.json')
	graphspace_utils.postGraph('classExampleEli','classExampleLab7.json','elispilio@reed.edu','shareabledrowssap')
	return
###anna's file parser for Distance matrices
def readFromFile(filename):
	"""
	parse an n-by-n matrix of distances and a header list.
	Input: name of file to parse
	Output: distance matrix D and header list
	"""
	D = []
	headers = []
	with open(filename) as fin:
		for line in fin:
			row = line.strip().split()
			if len(headers) == 0:
				headers = row
			else:
				row = [int(a) for a in row[1:]]
				D.append(row)
	return D,headers
##neighbor joining function

def neighb_join(D,headers):
	OTUs = len(D)+1
	int_points = 0
	T = []
	while int_points < OTUs -2:
		Q = []
		iterab = 0
		for pos1 in headers:
			for pos2 in headers:
				if pos1 != iterab:
					Q.append([])
					iterab = pos1
				if pos1 == pos2:
					Q[-1].append(100000)
				else:
					i = headers.index(pos1)
					j = headers.index(pos2)
					Q[-1].append((len(D)-2)*D[i][j] - sum(D[i]) - sum(D[j]))
		mini = (10000, [])
		for a in range(len(Q)):
			for b in range(len(Q)):
				if Q[a][b] < mini[0]:
					mini = (Q[a][b],[a,b])
		I = mini[1][0]
		J = mini[1][1]
		T.append([headers[I],headers[I] + headers[J],D[I][J]*.5 + (1/(2*len(D)-2))*(sum(D[I]) - sum(D[J]))])
		T.append([headers[J],headers[I]+headers[J],T[-1][2] - D[I][J]])
		newhead = headers[I]+headers[J]
		headers.append(headers[I] + headers[J])
		if I < J:
			del headers[J]
			del headers[I]
		else:
			del headers[I]
			del headers[J]
		col = []
		num = D[I][J]
		for row in D:
			new = .5*(row[I] + row[J] - num)
			row.append(new)
			col.append(new)
			if I < J:
				del row[J]
				del row[I]
			else:
				del row[I]
				del row[J]
		if I < J:
			del col[J]
			del col[I]
		else: 
			del col[I]
			del col[J]
		col.append(0)
		D.append(col)
		if J < I:
			del D[I]
			del D[J]
		else:
			del D[J]
			del D[I]
		int_points += 1
	return T

def getNodeAttributes(tree):
	nodes = []
	interior = []
	for a in tree:
		if a[0] not in nodes:
			nodes.append(a[0])
		if a[1] not in nodes:
			nodes.append(a[1])
	for b in nodes:
		if len(b) >1:
			interior.append(b)
	nodeAttrs = {}
	for n in nodes:
		nodeAttrs[n] = {}
		nodeAttrs[n]['id'] = n
		nodeAttrs[n]['content'] = n
		if n in interior:
			nodeAttrs[n]['background_color'] = '#ff0000'
		else:
			nodeAttrs[n]['background_color'] = '#0000ff'
	return nodeAttrs

def getEdgeAttributes(tree):
	attrs = {}
	edges = tree
	for e in edges:
		source = e[0]
		target = e[1]
		if source not in attrs:
			attrs[source] = {}
		attrs[source][target] = {}
		#attrs[source][target]['width'] = e[2]
		attrs[source][target]['popup'] = str(e[2])
	return attrs

if __name__ == '__main__':
	main()
