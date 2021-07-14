import pathlib
import os
import argparse
import picklify
import distance


def create_mer_count_files(organisms, k_mer_length, threads, dest_loc):
    # create the output folder if it does not exist
    dest_path = pathlib.Path(dest_loc)
    dest_path.mkdir(parents=True, exist_ok=True)

    for organism in organisms:
        if os.path.exists(f'{dest_loc}/{organism[0]}.txt'):
            if input(f"Reset {f'{dest_loc}/{organism[0]}.txt'}? [y/n]\n").upper() == 'Y':
                os.remove(f'{dest_loc}/{organism[0]}.txt')
            else:
                if input(f"Continue the program? [y/n]").upper() == 'N':
                    return

    for organism in organisms:
        print(f'Working on {organism[0]} files...')
        path = pathlib.Path(f"{organism[1]}/")

        for filename in path.iterdir():
            seqid = filename.stem
            os.system(f'jellyfish count -m {k_mer_length} -s 100M -t {threads} -C {filename}')
            os.system(f"jellyfish dump mer_counts.jf -c > dumped.txt")
            with open(f'{dest_loc}/{organism[0]}.txt', 'a') as k_mers_file:
                k_mers_file.write(f'>{seqid}\n')
                with open('dumped.txt', 'r') as dumped:
                    for line in dumped:
                        k_mers_file.write(f'{line}')
                    k_mers_file.write(f'\n')

    os.remove('dumped.txt')
    os.remove('mer_counts.jf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Count k-mers of k length from files and save them')
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('-k', type=int, help='the length of mers')
    mode.add_argument('-s', '--skip', dest='s', action='store_true', help='select this to skip kmerizing')
    parser.add_argument('-t', type=int, default=4, help='the amount of threads for jellyfish to use [default=4]')
    parser.add_argument('-o', '--out_dir', dest='o', default='./', help='location to save the files to [default=./]')
    parser.add_argument('-b', '--host_dir', dest='b', help='path to dir with bacteria (host) sequences')
    parser.add_argument('-v', '--vir_dir', dest='v', help='path to dir with virus sequences')
    parser.add_argument('-p', '--pickle', dest='p', action='store_true', help='select this to create pickle files')
    parser.add_argument('-d', '--distance', dest='d', 
        choices=['manhattan', 'cityblock', 'euclidean', 'canberra', 'chebyshev', 'cosine', 'braycurtis', 'pearson'],
                        help='choose which distance to calculate if any')
    args = parser.parse_args()

    if args.s:
        dest = args.o
    else:
        dest = f'{args.o}/k{args.k}'

    if not args.s:
        # dest = f'{args.o}/k{args.k}'
        organisms = []
        if args.v:
            organisms.append(['virus', args.v])
        if args.b:
            organisms.append(['host', args.b])
        create_mer_count_files(organisms, args.k, args.t, dest)

    if args.p:
        print(f'Creating pickle files...')
        # if args.s:
        #     dest = args.o
        # getting kmer_universe requires both host and virus files
        orgs = ['virus', 'host']

        # create the k_mer universe - a set of all kmers present across all the files for k_mers of this length
        kmer_universe = picklify.get_kmer_universe(dest, orgs)

        # kmer_dict - i stands for the index we're going to write the frequency of a particular kmer in an organism to
        kmer_dict = {kmer: i for i, kmer in enumerate(kmer_universe)}

        # create pickles
        picklify.create_pickles(dest, orgs, kmer_universe, kmer_dict)
    if args.d:
        print(f'Calculating {args.d} distance...')
        distance.distance(dest, args.d)

    print('DONE')
