import pickle
import heapq
from collections import defaultdict
from os.path import exists, join
from os import makedirs
from multiprocessing import Pool


class PPINetwork:
    def __init__(self, ppi_path, drug_path, output_dir, num_process):
        with open(ppi_path, "rb") as f:
            self.ppi = pickle.load(f)
        with open(drug_path, "rb") as f:
            self.drug_targets = pickle.load(f)
        self.drug_list = list(self.drug_targets.keys())
        if not exists(output_dir):
            makedirs(output_dir)
        self.output_dir = output_dir
        self.num_process = num_process
        print("total number of drugs: ", len(self.drug_list))

    def compute_distances_multi_drug_parallel(self, start_indx, end_index):
        with Pool(self.num_process) as p:
            p.map(self.compute_distances_single_drug, self.drug_list[start_indx:end_index + 1])

    def compute_distances_single_drug(self, src_drug):
        """
        Returns the shortest distances between src_drug's targets to the
        first-order neighbors of other drugs.
        """
        # all the target proteins of source drug
        targets = self.drug_targets[src_drug]
        targets = [x for x in targets if x.isdigit()]

        shortest_distances = {}
        for target in targets:
            # run dijkstra to find the shortest distances of all nodes to target
            min_dist = self.dijkstra(self.ppi, target)
            for dst_drug in self.drug_list:
                if dst_drug==src_drug:
                    continue
                dst_targets = self.drug_targets[dst_drug]
                dst_targets = [x for x in dst_targets if x.isdigit()]
                visited = set()
                for dst_target in dst_targets:
                    if dst_target in visited:
                        continue
                    visited.add(dst_target)
                    for _, neighbor in self.ppi[dst_target]:
                        path = src_drug + '-' + target + '-' + dst_drug + '-' + dst_target + '-' + neighbor
                        shortest_distances[path] = min_dist[neighbor]

        # save the results of src_drug
        out_file_path = join(self.output_dir, f"{src_drug}.pickle")
        with open(out_file_path, "wb") as f:
            pickle.dump(shortest_distances, f)

    @staticmethod
    def dijkstra(graph, start):
        """Returns the shortest distances of all other nodes in the graph."""
        visited = set()
    
        # shortest distance of each node to start
        min_dist = defaultdict(lambda: float('inf'))
        min_dist[start] = 0

        # use a heap to get the unvisited node that has shortest distance
        # heap element: (dist, node)
        heap = [(0, start)]

        # start dijkstra algorithm
        while heap:
            dist, node = heapq.heappop(heap)

            # Skip the stale shortest distances.
            if min_dist[node] < dist: continue

            visited.add(node)

            for weight, neighbor in graph[node]:
                if neighbor in visited: continue
                new_dist = dist + weight 

                # edge relaxation
                if new_dist < min_dist[neighbor]:
                    min_dist[neighbor] = new_dist
                    heapq.heappush(heap, (new_dist, neighbor))

        return min_dist


if __name__=="__main__":
    start_index, end_index = 70, 96
    ppi = PPINetwork(ppi_path="../ppi.pickle", drug_path="../drugs.pickle", output_dir="../output", num_process=2)
    ppi.compute_distances_multi_drug_parallel(start_index, end_index)