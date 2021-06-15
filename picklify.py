import numpy as np
import pickle
import argparse
import pathlib


def get_kmer_universe(path, organisms):
    kmer_universe = set()
    for organism in organisms:
        with open(f'./{str(path)}/{organism}.txt', 'r') as in_file:
            for line in in_file:
                if line.startswith('>'):
                    continue
                elif line.startswith('\n'):
                    continue
                else:
                    kmer_universe.add(line.split()[0])

    return kmer_universe


def create_pickles(path, organisms, kmer_universe, kmer_dict):
    # create pickle files with FREQUENCIES of kmers (guaranteed same order)
    for organism in organisms:
        # check if organism dir exists
        pathlib.Path(f'./{str(path)}/{organism}_pickles').mkdir(parents=True, exist_ok=True)

        with open(f'./{str(path)}/{organism}.txt', 'r') as in_file:
            current_org = {}
            curr_name = ''
            for line in in_file:
                if line.startswith('>'):
                    curr_name = line.strip()[1:]
                    current_org[curr_name] = {}
                elif line.startswith('\n'):
                    # create numpy array to write down to a pickle file
                    arr = np.zeros(len(kmer_universe))

                    # calculate the total amount of mers
                    total_mers = sum(current_org[curr_name].values())
                    # convert mer counts to mer frequencies
                    for mer in current_org[curr_name]:
                        current_org[curr_name][mer] = current_org[curr_name][mer] / total_mers
                        idx = kmer_dict[mer]
                        arr[idx] = current_org[curr_name][mer]

                    # save the array to a pickle file
                    with open(f'./{str(path)}/{organism}_pickles/{curr_name}.p', 'wb') as oh:
                        pickle.dump(arr, oh)

                    # reset the current org values
                    current_org = {}
                    curr_name = ''
                else:
                    segs = line.strip().split()
                    current_org[curr_name][segs[0]] = int(segs[1])


def calculate_manhattan_distance(path):
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
    parser = argparse.ArgumentParser('Compare the k-mer frequencies between phages and hosts')
    parser.add_argument('in_dir', help='path to dir with pre-generated input files')
    args = parser.parse_args()

    path = pathlib.Path(f'./{args.in_dir}/')

    organisms = ['virus', 'host']

    # create the k_mer universe - a set of all kmers present across all the files for k_mers of this length
    kmer_universe = get_kmer_universe(path, organisms)

    # kmer_dict - i stands for the index we're going to write the frequency of a particular kmer in an organism to
    kmer_dict = {kmer: i for i, kmer in enumerate(kmer_universe)}

    # create pickles
    create_pickles(path, organisms, kmer_universe, kmer_dict)

