import pathlib
import argparse
import numpy as np
import pickle


def manhattan(path):
    # read pickles and create the results file
    virus_path = pathlib.Path(f'./{str(path)}/virus_pickles/')
    host_path = pathlib.Path(f'./{str(path)}/host_pickles/')
    with open(f'{str(path)}/results.txt', 'w') as results_file:
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
                    results_file.write(f'{virus[0]}\t{host_name}\t'
                                       f'{np.sum(np.abs(virus[1] - h_list)) / h_list.size}\n')


if __name__ == '__main__':
    # enable arguments from terminal
    parser = argparse.ArgumentParser('Compare the k-mer frequencies between phages and hosts')
    parser.add_argument('in_dir', help='path to dir with pre-generated input files')
    parser.add_argument('distance', choices=['manhattan'],
                        help='choose which distance to calculate if any')
    args = parser.parse_args()
    path = pathlib.Path(f'./{args.in_dir}/')

    # TODO - add possibility to choose which kind of organism to use
    # TODO - add different types of distance

    print(f'Calculating {args.distance} distance...')
    manhattan(path)
