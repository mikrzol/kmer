import pathlib
import argparse
import numpy as np
import pickle
import scipy.spatial as sp

'''
# tutaj funkcje: manhattan, euclidean itp. - zwracają konkretną wartość float
def manhattan(virus, h_list):
    return np.sum(np.abs(virus - h_list))


def euclidean(virus, h_list):
    # euclidean distance is the i2 norm, so we can use the numpy.linalg.norm function
    return np.linalg.norm(virus - h_list)


def canberra(virus, h_list):
    return sp.distance.canberra(virus, h_list)


def chebyshev(virus, h_list):
    return sp.distance.chebyshev(virus, h_list)


def cosine(virus, h_list):
    return sp.distance.cosine(virus, h_list)
'''

def distance(path, dist_type):
    # read pickles and create the results file
    virus_path = pathlib.Path(f'./{str(path)}/virus_pickles/')
    host_path = pathlib.Path(f'./{str(path)}/host_pickles/')
    with open(f'{str(path)}/{dist_type}_results.txt', 'w') as results_file:
        # load all the virus pickles to memory
        virus_list = []
        for virus_file in virus_path.iterdir():
            virus_name = virus_file.stem
            with open(virus_file, 'rb') as vh:
                virus_list.append((virus_name, pickle.load(vh)))

        for host_file in host_path.iterdir():
            host_name = host_file.stem
            with open(host_file, 'rb') as hh:
                h_list = pickle.load(hh)
                # TODO maybe - add multiprocessing - calculate distances,
                #  save them to a list and write the whole list to a file
                for virus in virus_list:
                    # tutaj można skorzystać z yield
                    # virus[0] -> name; virus[1] -> pickled array 
                    results_file.write(f'{virus[0]}\t{host_name}\t'
                                        f'{eval(f"sp.distance.{dist_type}(virus[1], h_list)")}\n'
                                       # f'{eval(dist_type + f"(virus[1], h_list)")}\n'
                                       )


if __name__ == '__main__':
    # enable arguments from terminal
    parser = argparse.ArgumentParser(
        'Compare the k-mer frequencies between phages and hosts')
    parser.add_argument(
        'in_dir', help='path to dir with pre-generated input files')
    parser.add_argument('distance', 
        choices=['manhattan', 'euclidean', 'canberra', 'chebyshev', 'cosine', 'braycurtis'],
                        help='choose which distance to calculate if any')
    args = parser.parse_args()
    path = pathlib.Path(f'./{args.in_dir}/')
    args.distance = 'cityblock' if args.distance == 'manhattan' else args.distance

    # TODO - add possibility to choose which kind of organism to use
    # TODO - add different types of distance

    print(f'Calculating {args.distance} distance...')
    distance(path, args.distance)

    # TODO współczynnik korelacji liniowej - Pearsona - trzeba sprytnie zmienić na dystans
    # wyżej: najmniejszy dystans powinna mieć 1, a najmniejszy: -1
