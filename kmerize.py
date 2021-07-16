import pathlib
import os
import argparse
import picklify
import distance


def create_mer_count_files(organisms, k_mer_length, threads, dest_loc):
    print(f'dest_loc = {dest_loc}')
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
            os.system(f'jellyfish count -m {k_mer_length} -s 100M -t {threads} -C {filename} -o {dest_loc}/mer_counts.jf')
            os.system(f"jellyfish dump {dest_loc}/mer_counts.jf -c > {dest_loc}/dumped.txt")
            with open(f'{dest_loc}/{organism[0]}.txt', 'a') as k_mers_file:
                k_mers_file.write(f'>{seqid}\n')
                with open(f'{dest_loc}/dumped.txt', 'r') as dumped:
                    for line in dumped:
                        k_mers_file.write(f'{line}')
                    k_mers_file.write(f'\n')

    os.remove(f'{dest_loc}/dumped.txt')
    os.remove(f'{dest_loc}/mer_counts.jf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Count k-mers of k length from files and save them')
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('-k', type=int, help='the length of mers')
    mode.add_argument('-s', '--skip', dest='s', action='store_true', help='select this to skip kmerizing')
    parser.add_argument('-t', type=int, default=4, help='the amount of threads for jellyfish to use [default=4]')
    parser.add_argument('-o', '--out_dir', dest='o', default='./', help='location to save the files to [default=./]')
    parser.add_argument('-b', '--host_dir', dest='b', help='path to dir with bacteria (host) sequences')
    parser.add_argument('-v', '--vir_dir', dest='v', help='path to dir with virus sequences')
    args = parser.parse_args()

    dest = args.o

    if not args.s:
        organisms = []
        if args.v:
            organisms.append(['virus', args.v])
        if args.b:
            organisms.append(['host', args.b])
        create_mer_count_files(organisms, args.k, args.t, dest)

    print('DONE')
