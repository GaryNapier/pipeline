from ete3 import Tree, ClusterTree, TreeStyle, NodeStyle, AttrFace, ProfileFace, TextFace, faces
from ete3.treeview.faces import add_face_to_node

# ------------------------------
# Tree and clustering functions
# ------------------------------

class tree:
    def __init__(self,treefile):
        self.treefile = treefile
        self.t = ete3.Tree(treefile)

    def distance_from_root(self):
        dists = {}
        for leaf in tqdm(self.t.get_leaves()):
            dists[leaf.name] = self.t.get_distance(leaf)
        return dists

    def extract_distance_matrix(self):
        return extract_distance_matrix(self)

def read_tree(tree_file):
	with open(tree_file, mode='r') as newick_file:
		t = Tree(newick_file.read())
	# Midpoint root the tree
	# Calculate the midpoint node
	md_pt = t.get_midpoint_outgroup()
	# and set it as tree outgroup
	t.set_outgroup(md_pt)
	# Write midpoint rooted as newick file
	with open(tree_file+".rooted", mode='w') as newick_file:
		newick_file.write(t.write())
	return t

def get_clusters(dists, prefix, cutoff=10, remove_singletons=False):
	# dists = self.get_plink_dist()
	edges = []
	tmp_node_set = set()
	samples = list(dists.columns)
	transmission_samples = []
	for row_ind, row in enumerate(dists.index.values):
		for col_ind, col in enumerate(dists.columns.values):
			if col >= row: # Lower triangle
				continue
			#if subdist_within_table.loc[row, col] <= outliers_within_stat[0] or subdist_within_table.loc[row, col] >= outliers_within_stat[1]:
			if dists.loc[row, col] <= cutoff:
				edge = {"source": samples[row_ind], "target": samples[col_ind], "snps": dists.loc[row, col]}
				tmp_node_set.add(samples[row_ind])
				tmp_node_set.add(samples[row_ind])
				edges.append(edge)

	nodes = [{"id": s} for s in tmp_node_set] if remove_singletons else [
		{"id": s} for s in samples]
	graph = {"nodes": nodes, "edges": edges}
	return graph

class transmission_graph:
	# def __init__(self, filename):
	def __init__(self, graph):
		# tmp = json.load(open(filename))
		# self.json_graph = tmp
		self.json_graph = graph
		self.graph = networkx.Graph()
		# self.graph.add_nodes_from([x["id"] for x in tmp["nodes"]])
		self.graph.add_nodes_from([x["id"] for x in self.json_graph["nodes"]])
		for edge in self.json_graph["edges"]:
			self.graph.add_edges_from([(edge["source"], edge["target"])])
		self.clusters = sorted(list(connected_components(self.graph)), key=lambda x: len(x), reverse=True)

	def add_meta_data(self, csvfile):
		for row in csv.DictReader(open(csvfile)):
			found = False
			for i, node in enumerate(self.json_graph["nodes"]):
				if node["id"] == row["id"]:
					found = True
					break
			if found:
				for column in set(row.keys())-set(["id"]):
					if column not in self.json_graph["nodes"][i]:
						self.json_graph["nodes"][i][column] = row[column]

	def extract_clusters(self):
		print(self.clusters)

def get_leaves_in_clusts(tree_file, node_cutoff, bootstrap_cutoff):

    # Read tree file
    t = read_tree(tree_file)
    # Parse tree and get lists of samples for each clade ('leaves')
    # Clades have to be > 50 bootstrap an min of 20 samples (by default)
    leaves = list()
    for node in t.traverse("preorder"):
        # Do some analysis on node
        if len(node) > node_cutoff and node.support >= bootstrap_cutoff:
            # print(node)
            # print(node.name)
            # print(node.support)
            leaves.append([n.name for n in node.get_leaves()])
    return leaves

def extract_distance_matrix(t):
    dists = {}
    for i,leafi in enumerate(t.get_leaves()):
        dists[leafi.name] = {}
        for j,leafj in tqdm(enumerate(t.get_leaves())):
            if i>=j: continue
            dists[leafi.name][leafj.name] = leafi.get_distance(leafj)
    return dists

def main_extract_leaf_root_dists(args):
    t = tree(args.tree)
    dists = t.distance_from_root()
    for sample in dists:
        print(f"{sample}\t{dists[sample]}")

def main_extract_paired_dists(args):
    t = tree(args.tree)
    print(t.extract_distance_matrix())

def main_subsample(args):
    t = tree(args.tree)
    dists = t.distance_from_root()
    median_dist = median(list(dists.values()))
    dist_cutoff = args.arg1*median_dist
    for node in t.t.traverse():
        farthest_dist = node.get_farthest_leaf()[1]
        if farthest_dist > dist_cutoff: continue
        dists = extract_distance_matrix(node)
        print(dists)