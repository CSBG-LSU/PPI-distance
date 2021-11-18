from ppi_dijkstra import PPINetwork
from os import listdir, makedirs
from os.path import isfile, join, exists
import pickle
from re import split
from collections import defaultdict
from random import choices
from multiprocessing import Pool


class PPINetworkRandom(PPINetwork):
    def __init__(self, ppi_path, drug_path, output_dir, output_db_dir, random_output_dir, num_process):
        super().__init__(ppi_path, drug_path, output_dir, output_db_dir, num_process)
        if not exists(random_output_dir):
            makedirs(random_output_dir)
        self.random_output_dir = random_output_dir

        # compute number of random destinations for each node
        self.nodes = list(self.ppi.keys())
        #num_node_pairs = self.compute_num_node_pairs()
        #self.num_random_per_node = int(num_node_pairs / len(self.nodes))

        # pre-computed using self.compute_num_node_pairs
        self.num_random_per_node = 4131489
        
        #print(f"number of nodes: {len(self.nodes)}")
        #print(f"number of random pairs: {num_node_pairs}")
        #print(f"number of random destination nodes per node: {self.num_random_per_node}")

    def compute_random_distances_parallel(self):
        with Pool(self.num_process) as p:
            p.map(self.compute_random_distances_single_node, self.nodes)
        
    def compute_random_distances_single_node(self, node):
        """Compute distances from input node to random nodes on the graph"""
        random_nodes = choices(self.nodes, k=self.num_random_per_node)
        min_dist = self.dijkstra(self.ppi, node)
        result = [min_dist[random_node] for random_node in random_nodes]
        with open(join(self.random_output_dir, f"{node}.pickle"), "wb") as f:
            pickle.dump(result, f)

    def compute_num_node_pairs(self):
        num_node_pairs = defaultdict(lambda: 0)
        out_files = [f for f in listdir(self.output_dir) if isfile(join(self.output_dir, f))]
        # count the number of items in each output file
        for out_file in out_files:
            node = split('-|\.', out_file)[1]
            with open(join(self.output_dir, out_file), "rb") as f:
                paths = pickle.load(f)
            num_node_pairs[node] += len(paths)
        return sum(num_node_pairs.values())


if __name__=="__main__":
    ppi = PPINetworkRandom(
        ppi_path="../ppi.pickle", 
        drug_path="../drugs.pickle", 
        output_dir="../output/", 
        output_db_dir="../output_db/",
        random_output_dir="../output_random/",
        num_process=28
    )
    ppi.compute_random_distances_parallel()