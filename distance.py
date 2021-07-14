import pathlib
import argparse
import numpy as np
import pickle
import math
import time
# import scipy.spatial as sp


# TODO implement these functions myself
def manhattan(virus, host):
    return np.sum(np.abs(virus - host))


def euclidean(virus, host):
    # euclidean distance is the i2 norm, so we can use the numpy.linalg.norm function
    return np.linalg.norm(virus - host)


def chebyshev(virus, host):
    return np.max(np.abs(virus - host))


def canberra(virus, host):
    return np.nansum(abs(virus - host) / (abs(virus) + abs(host)))


def cosine(virus, host):
    # it's the distance, so we need to subtract the similarity coeff from 1
    a = np.dot(virus, host)
    b = np.sqrt(np.dot(virus, virus)) * np.sqrt(np.dot(host, host))
    return 1 - np.divide(a, b, out=np.zeros_like(a), where=b!=0)


def braycurtis(virus, host):
    a = np.sum(np.abs(virus - host))
    b = np.sum(np.abs(virus + host))
    return np.divide(a, b, out=np.zeros_like(a), where=b!=0)


def pearson(virus, host):
    # return 1 - np.corrcoef(virus, host)[0][1]
    v_avg = np.sum(virus) / len(virus)
    h_avg = np.sum(host) / len(host)
    numerator = np.sum(np.subtract(virus, v_avg) * np.subtract(host, h_avg))
    denominator = np.sqrt(np.sum((virus - v_avg)**2)) * np.sqrt(np.sum((host - h_avg)**2))
    # need to return 1 - coeff -> the distance
    return 1 - np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator!=0)


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
                                        f'{eval(dist_type + f"(virus[1], h_list)")}\n'
                                       # f'{eval(f"sp.distance.{dist_type}(virus[1], h_list)")}\n'
                                       )


if __name__ == '__main__':
    # enable arguments from terminal
    parser = argparse.ArgumentParser(
        'Compare the k-mer frequencies between phages and hosts')
    parser.add_argument(
        'in_dir', help='path to dir with pre-generated input files')
    parser.add_argument('distance', 
        choices=['manhattan', 'euclidean', 'canberra', 'chebyshev', 'cosine', 'braycurtis', 'pearson'],
                        help='choose which distance to calculate if any')
    args = parser.parse_args()
    path = pathlib.Path(f'./{args.in_dir}/')
    args.distance = 'cityblock' if args.distance == 'manhattan' else args.distance

    # TODO - add possibility to choose which kind of organism to use
    # TODO - add different types of distance

    print(f'Calculating {args.distance} distance...')
    start = time.time()
    distance(path, args.distance)
    end = time.time()
    print('DONE')
    print(f'Time elapsed: {end - start:.2f}')

    # TODO współczynnik korelacji liniowej - Pearsona - trzeba sprytnie zmienić na dystans
    # wyżej: najmniejszy dystans powinna mieć 1, a najmniejszy: -1
