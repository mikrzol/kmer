import pathlib
import argparse
import numpy as np
from pickle import load as p_load
import time


# TODO implement these functions myself
def cityblock(virus, host):
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
    # custom implementation is faster than numpy's
    v_avg = np.sum(virus) / len(virus)
    h_avg = np.sum(host) / len(host)
    numerator = np.sum(np.subtract(virus, v_avg) * np.subtract(host, h_avg))
    denominator = np.sqrt(np.sum((virus - v_avg)**2)) * np.sqrt(np.sum((host - h_avg)**2))
    # need to return 1 - coeff -> the distance
    return 1 - np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator!=0)


def distance(path, dist_type):
    dist_type = 'cityblock' if dist_type == 'manhattan' else dist_type
    # read pickles and create the results file
    virus_path = pathlib.Path(f'./{str(path)}/virus_pickles/')
    host_path = pathlib.Path(f'./{str(path)}/host_pickles/')
    with open(f'{str(path)}/{dist_type}_results.txt', 'w') as results_file:
        # load all the virus pickles to memory
        virus_list = []
        for virus_file in virus_path.iterdir():
            virus_name = virus_file.stem
            with open(virus_file, 'rb') as vh:
                virus_list.append((virus_name, p_load(vh)))

        for host_file in host_path.iterdir():
            host_name = host_file.stem
            with open(host_file, 'rb') as hh:
                h_list = p_load(hh)
                for virus in virus_list:
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
        choices=['manhattan', 'cityblock', 'euclidean', 'canberra', 'chebyshev', 'cosine', 'braycurtis', 'pearson'],
                        help='choose which distance to calculate if any')
    args = parser.parse_args()
    path = pathlib.Path(f'./{args.in_dir}/')

    print(f'Calculating {args.distance} distance...')
    start = time.time()
    distance(path, args.distance)
    end = time.time()
    print('DONE')
    print(f'Time elapsed: {end - start:.2f}s')
