import argparse
import kmerize
import distance
import picklify

if __name__ == '__main__':
    parser = argparse.ArgumentParser('mikrz kmerize control panel')
    parser.add_argument('-k', '--kmerize', type=int,
        help='pass in the desired length of mers to perform kmerize step (create mer count files)')
    parser.add_argument('-o', '--out_dir', dest='o', default='./', help='location to save the files to [default=./]'),
    parser.add_argument('-b', '--host_dir', dest='b', help='path to dir with bacteria (host) sequences')
    parser.add_argument('-v', '--vir_dir', dest='v', help='path to dir with virus sequences')
    parser.add_argument('-t', type=int, default=4, help='the amount of threads for jellyfish to use (for kmerize step) [default=4]')
    parser.add_argument('-p', '--picklify', 
        action='store_true', help='select this to picklify (create pickle files)')
    parser.add_argument('-d', '--distance', 
        choices=['manhattan', 'cityblock', 'euclidean', 'canberra', 'chebyshev', 'cosine', 'braycurtis', 'pearson'],
        help='choose which distance to calculate if any')
    args = parser.parse_args()

    if args.kmerize:
        organisms = []
        if args.v:
                organisms.append(['virus', args.v])
        if args.b:
            organisms.append(['host', args.b])
        kmerize.create_mer_count_files(organisms, args.kmerize, args.t, args.o)
    
    if args.picklify:
        print(f'Creating pickle files...')
        # getting kmer_universe requires both host and virus files
        orgs = ['virus', 'host']

        # create the k_mer universe - a set of all kmers present across all the files for k_mers of this length
        kmer_universe = picklify.get_kmer_universe(args.o, orgs)

        # kmer_dict - i stands for the index we're going to write the frequency of a particular kmer in an organism to
        kmer_dict = {kmer: i for i, kmer in enumerate(kmer_universe)}

        # create pickles
        picklify.create_pickles(args.o, orgs, kmer_universe, kmer_dict)
   
    if args.distance:
        print(f'Calculating {args.distance} distance...')
        distance.distance(args.o, args.distance)
    
    print('DONE')
