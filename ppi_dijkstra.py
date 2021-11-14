import pickle
import heapq
from collections import defaultdict
from os.path import exists, join
from os import makedirs
from multiprocessing import Pool
import sqlite3


class PPINetwork:
    def __init__(self, ppi_path, drug_path, output_dir, output_db_dir, num_process):
        with open(ppi_path, "rb") as f:
            self.ppi = pickle.load(f)
        with open(drug_path, "rb") as f:
            self.drug_targets = pickle.load(f)
        self.drug_list = list(self.drug_targets.keys())
        if not exists(output_dir):
            makedirs(output_dir)
        self.output_dir = output_dir
        if not exists(output_db_dir):
            makedirs(output_db_dir)
        self.output_db_dir = output_db_dir
        self.num_process = num_process
        print("total number of drugs: ", len(self.drug_list))

    def compute_distances_multi_drug_parallel_db(self, start_indx, end_index):
        with Pool(self.num_process) as p:
            p.map(self.compute_distances_single_drug_db, self.drug_list[start_indx:end_index + 1])

    def compute_distances_single_drug_db(self, src_drug):
        """
        Saves the shortest distances between src_drug's targets to the
        first-order neighbors of other drugs to a database.
        """
        # initialize a database for each drug
        database_path = join(self.output_db_dir, f"{src_drug}.db")
        con, cur = self.init_db(database_path)

        # all the target proteins of source drug
        targets = self.drug_targets[src_drug]
        targets = [x for x in targets if x.isdigit()]

        for target in targets:
            # run dijkstra to find the shortest distances of all nodes to target
            min_dist = self.dijkstra(self.ppi, target)
            for dst_drug in self.drug_list:
                if dst_drug==src_drug:
                    continue
                dst_targets = self.drug_targets[dst_drug]
                dst_targets = [x for x in dst_targets if x.isdigit()]
                for dst_target in dst_targets:
                    for _, neighbor in self.ppi[dst_target]:
                        path = src_drug + '-' + target + '-' + dst_drug + '-' + dst_target + '-' + neighbor
                        cur.execute(f"INSERT INTO shortest_distances VALUES ('{path}', {min_dist[neighbor]})")
                        con.commit()

        # remove duplicate rows
        cur.execute('''DROP TABLE IF EXISTS temp_table''')
        cur.execute('''CREATE TABLE temp_table as SELECT DISTINCT * FROM shortest_distances''')
        cur.execute('''DELETE FROM shortest_distances''')
        cur.execute('''INSERT INTO shortest_distances SELECT * FROM temp_table''')
        cur.execute('''DROP TABLE IF EXISTS temp_table''')
        con.commit()
        con.close()

    @staticmethod
    def init_db(db_path):
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur.execute('''DROP TABLE IF EXISTS shortest_distances''')
        cur.execute('''CREATE TABLE shortest_distances
               (path text, distance real)''')
        con.commit()
        return con, cur

    def compute_distances_multi_drug_parallel(self, start_indx, end_index):
        with Pool(self.num_process) as p:
            p.map(self.compute_distances_single_drug, self.drug_list[start_indx:end_index + 1])

    def compute_distances_single_drug(self, src_drug):
        """
        Saves the shortest distances between src_drug's targets to the
        first-order neighbors of other drugs to pickle files.
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
                for dst_target in dst_targets:
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
    #start_index, end_index = 70, 96
    #start_index, end_index = 50, 69
    #start_index, end_index = 20, 49
    start_index, end_index = 0, 96
    ppi = PPINetwork(
        ppi_path="../ppi.pickle", 
        drug_path="../drugs.pickle", 
        output_dir="../output", 
        output_db_dir="../output_db/",
        num_process=16
    )
    #ppi.compute_distances_multi_drug_parallel(start_index, end_index)
    ppi.compute_distances_multi_drug_parallel_db(start_index, end_index)