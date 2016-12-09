import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import itertools
G = nx.dense_gnm_random_graph(27,42)
nx.enumerate_all_cliques(G)
degrees = []
for n in G.nodes():
	degree = len(G.neighbors(n))
	G.node[n]['degree'] = degree
	degrees.append(degree)
pos = nx.spring_layout(G)
normalizer = (max(degrees)-min(degrees))
def rgb_to_hex(num):
	r = int(num*255)
	b = int((1-num)*255)
	g = int(0*255)
	return '#{:02x}{:02x}{:02x}'.format(r,g,b)
for node in G.nodes():
	nx.draw_networkx_nodes(G,pos,nodelist=[node],node_color=rgb_to_hex((G.node[node]['degree']-min(degrees)) / normalizer),node_size = 500,alpha=0.8)
nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
for sub in list(nx.find_cliques(G)):
	if len(sub) > 2:
		for pair in itertools.combinations(sub,2):
			newpair=[pair[0],pair[1]]
			print(newpair)
			nx.draw_networkx_edges(G,pos,edgelist=[newpair],width=1.0,alpha=0.5,edge_color="#FF0000")
plt.savefig("CliquesAndDegreeColor.png")