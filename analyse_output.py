from os import listdir
from os.path import isfile, join
import pickle


def get_output_list(output_dir):
    output_files = [f for f in listdir(output_dir) if isfile(join(output_dir, f))]
    return output_files

def read_output_file(output_dir, src_drug, src_target):
    file_name = f"{src_drug}-{src_target}.pickle"
    file_path = join(output_dir, file_name)
    with open(file_path, "rb") as f:
        dist_dict = pickle.load(f)
    return dist_dict

def query_distance(dist_dict, src_drug, src_target, dst_drug, dst_target, dst_target_neighbor):
    """Read the shortest distance to a first-order neighbor of dst_target."""
    path = f"{src_drug}-{src_target}-{dst_drug}-{dst_target}-{dst_target_neighbor}"
    return dist_dict[path]


if __name__ == "__main__":
    pass