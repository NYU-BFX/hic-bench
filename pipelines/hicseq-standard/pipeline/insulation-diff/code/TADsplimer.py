import sys
import argparse
import os
import numpy as np
import similarity_score
import TAD_calling
import TAD_split


def printHelp():
    print('\nTADsplimer version 1.0.2')
    print('For help information for each function, try:\npython3 TADsplimer.py <function> -h')
    print('\nFunctions:')
    print('\tsplit_TADs:\n\t\tsplit TAD detection(input contact maps in two conditions)\n')
    print('\tTAD_calculator:\n\t\tidentify the topological domain\n')
    print('\tsplit_TADs_alternate:\n\t\tsplit TAD detection (input contact maps in two conditions with pre-called TADs)\n')
    print('\tTAD_similarity:\n\t\tfour algorithms of TAD similarity\n')
    print('')


def TAD_calculator(command='TAD_calculator'):
    '''
    Description:
        This function provide an entrance to identify TAD.
    parameters:
        none
    '''

    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        print("\nusage:\npython3 TADsplimer.py TAD_calculator <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python3 TADsplimer.py TAD_calculator -h\n")
        return 0

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython3 TADsplimer.py TAD_calculator "
                                           "[optional arguments]\n\n", description='')
    parser.add_argument('command', default=None,
                        help="set as 'TAD_calculator' to identify TADs relative to Hi-C contact map")

    parser.add_argument('-c', '--contact_map', dest="contact_map", default=None,
                        help="path to Hi-C contact map")

    parser.add_argument('-u', '--up', dest="up", default=0, type=float,
                        help="up cutoff for Hi-C contact map")

    parser.add_argument('-d', '--down', dest="down", default=0, type=float,
                        help="down cutoff for Hi-C contact map")

    parser.add_argument('-o', '--TAD_output', dest="output", default=None,
                        help="path for the output file of TADs")

    parser.add_argument('-p', '--TAD_plot', dest="plot", default=1, type=int,
                        help="Set to 1 to plot the contact map and TAD, else set to 0 to cancel this analysis.")

    parser.add_argument('--sub_map', dest="sub_map", default=1, type=int,
                        help="Set to 1 to output sub contact maps and TADs, else set to 0 to cancel this analysis.")

    args = parser.parse_args()

    if not os.path.exists(args.output):
        print('path does not exist\n')
    else:
        up, down, _ = TAD_calling.file_split(args.contact_map, args.output)
        dirs = os.listdir(args.output)
        if args.up != 0 or args.down != 0:
            up = args.up
            down = args.down
        tad = np.empty(shape=[0, 3])
        for file in dirs:
            if file[-7:] == '.subchr':
                start = file.split('.')[-3]
                tad_temp = TAD_calling.TAD_calling(os.path.join(args.output, file), up, down)
                tad_location = tad_temp
                if args.plot == 1:
                    try:
                        TAD_calling.TAD_plot(os.path.join(args.output, file), os.path.join(args.output, file + '.tiff'),
                                             tad_temp, down, up)
                    except:
                        pass

                try:
                    tad_location[:, 1:3] = tad_location[:, 1:3] + int(start)
                    if not np.array_equal(tad_location, np.array([0, int(start), int(start)])):
                        tad = np.vstack((tad, tad_location))
                    if args.sub_map == 1:
                        np.savetxt(os.path.join(args.output, file + '.TAD.txt'), tad_location, delimiter="\t")
                except:
                    tad_location = np.array([[0, 0, 0]])
                    if args.sub_map == 1:
                        np.savetxt(os.path.join(args.output, file + '.TAD.txt'), tad_location, delimiter="\t")
        sorted_idx = np.lexsort(tad[:, 1:3].T)
        sorted_tad = tad[sorted_idx,:]
        row_mask = np.append([True], np.any(np.diff(sorted_tad[:, 1:3], axis=0), 1))
        unique_tad = sorted_tad[row_mask]
        np.savetxt(os.path.join(args.output, file + '.TAD.merge.txt'), unique_tad, delimiter="\t")


def split_TADs_alternate(command='split_TADs_alternate'):
    '''
    corner split algorithm for identifying split TAD
    '''

    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython3 TADsplimer.py TAD_calculator <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python3 TADsplimer.py TAD_calculator -h\n")
        return 0

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython3 TADsplimer.py corner_split <contact_map_file_paths> "
                                           "[optional arguments]\n\n", description='')
    parser.add_argument('command', default=None,
                        help="set as 'split_TADs_alternate' to identify split TAD")

    parser.add_argument('-c', '--contact_maps', dest="contact_map", default=None,
                        help="paths to Hi-C contact maps in two conditions. paths must be separated by the comma ','.")

    parser.add_argument('--contact_maps_aliases', dest="aliases", default=None,
                        help="A set of short aliases for the contact map. Paths must be separated by the comma ','.")

    parser.add_argument('-t', '--TAD', dest="TAD", default=None,
                        help="input files of TADs for two compared Hi-C contact maps. Paths must be separated by the"
                             " comma ','.")

    parser.add_argument('-u', '--up_cutoff', dest="up", default="0,0",
                        help="up cutoff for two compared Hi-C contact maps, paths must be separated by the comma ','.")

    parser.add_argument('-j', '--adjust_quality', dest="adjust_quality", default=0, type=int,
                        help="set as 1 to normalize sequence quality for two Hi-C contact maps, set as 0 not to "
                             "normalize sequence quality for two Hi-C contact maps")

    parser.add_argument('-o', '--output', dest="output", default=None,
                        help="path to output files")

    parser.add_argument('-d', '--split_direction', dest="direction", default=0, type=int,
                        help="set as 0: output TADs split in both two contact maps, set as 1: output TADs split in "
                             "contact map 1, set as 2: output TADs split in contact map 2")

    args = parser.parse_args()

    file = args.contact_map.split(',')
    TAD = args.TAD.split(',')
    up = args.up.split(',')
    aliases = args.aliases.split(',')

    TAD1 = np.loadtxt(TAD[0])
    TAD2 = np.loadtxt(TAD[1])

    TAD_s, TAD1_only, TAD2_only= TAD_split.TAD_matching(TAD2, TAD1)

    D1, D2, D3 = TAD_split.corner_split_score(TAD1, TAD2, TAD_s)

    map1 = np.loadtxt(file[0])
    map2 = np.loadtxt(file[1])

    up1, _ , _ = TAD_calling.file_split(TAD[0], outfile="", print_subcontact=1)
    up2, _ , _ = TAD_calling.file_split(TAD[1], outfile="", print_subcontact=1)

    if up[0] != 0 or up[1] != 0:
        up1 = float(up[0])
        up2 = float(up[1])

    if args.adjust_quality == 0:
        split1_1, split2_1, loc_d, loc_u = TAD_split.split_region(file[1], file[0], D3, D1, up2, up1, 1)
    elif args.adjust_quality == 1:
        map1 = np.loadtxt(file[0])
        map2 = np.loadtxt(file[1])
        ratio1 = TAD_split.get_ratio(map1, TAD1, float(up[0]))
        ratio2 = TAD_split.get_ratio(map2, TAD2, float(up[1]))
        fold = np.mean(ratio1)/np.mean(ratio2)
        split1_1, split2_1, loc_d, loc_u = TAD_split.split_region(file[1], file[0], D3, D1, up2, up1, fold)
    else:
        print("\nusage:\npython3 TADsplimer.py TAD_calculator <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python3 TADsplimer.py TAD_calculator -h\n")
    try:
        if split1_1 == 0:
            print('No split TAD')
    except:
        if (not os.path.exists(args.output)):
            os.makedirs(args.output)
        else:
            scc = similarity_score.similarity_scc(map1, map2, loc_u[:, 0:3])
            Laplacian = similarity_score.similarity_Laplacian(map1, map2, loc_u[:, 0:3])
            hash = similarity_score.hash_similarity(map1, map2, loc_u[:, 0:3])
            loc_u2 = np.c_[loc_u, scc[:, 3], Laplacian[:, 3], hash[:, 3]]
            np.savetxt(os.path.join(args.output, aliases[0]+'->'+aliases[1]+'.split.txt'), loc_d)
            np.savetxt(os.path.join(args.output, aliases[0]+'->'+aliases[1]+'.merge.txt'), loc_u2)


def split_TADs(command='split_TADs'):
    '''
    corner split algorithm for identifying split TAD
    '''

    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython3 TADsplimer.py TAD_calculator <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python3 TADsplimer.py TAD_calculator -h\n")
        return 0

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython3 TADsplimer.py corner_split <contact_map_file_paths> "
                                           "[optional arguments]\n\n", description='')
    parser.add_argument('command', default=None,
                        help="set as 'split_TADs' to identify split TADs")

    parser.add_argument('-c', '--contact_maps', dest="contact_map", default=None,
                        help="paths to two contact maps. paths must be separated by the comma ','.")

    parser.add_argument('--contact_maps_aliases', dest="aliases", default=None,
                        help="A set of short aliases for two contact maps. Paths must be separated by the comma ','.")

    parser.add_argument('--up_cutoff', dest="up_cutoff", default="None",
                        help="paths for up cutoff of two compared Hi-C contact maps, paths must be separated by the comma ','.")

    parser.add_argument('--down_cutoff', dest="down_cutoff", default="None",
                        help="paths for down cutoff of two compared Hi-C contact maps, paths must be separated by the comma ','.")

    parser.add_argument('-j', '--adjust_quality', dest="adjust_quality", default=0, type=int,
                        help="set as 1 to normalize sequence quality for two Hi-C contact maps, set as 0 not to "
                             "normalize sequence quality for two Hi-C contact maps")

    parser.add_argument('-o', '--output', dest="output", default=None,
                        help="path to output files")

    parser.add_argument('-d', '--split_direction', dest="direction", default=0, type=int,
                        help="set as 0: output TADs split in both two contact maps, set as 1: output TADs split in "
                             "contact map1, set as 2: output TADs split in contact map2")

    args = parser.parse_args()

    file = args.contact_map.split(',')
    aliases = args.aliases.split(',')

    tmp_dir1 = os.path.join(args.output, aliases[0])
    if not os.path.isdir(tmp_dir1):
        os.makedirs(tmp_dir1)

    up1, down1, files1 = TAD_calling.file_split(file[0], outfile=tmp_dir1, print_subcontact=0, aliases=aliases[0])

    print('auto:', up1, down1)

    if args.up_cutoff != 'None' and args.down_cutoff != 'None':
        up_cutoff = args.up_cutoff.split(',')
        up1 = float(up_cutoff[0])
        down_cutoff = args.down_cutoff.split(',')
        down1 = float(down_cutoff[0])
        print('manual:', up1, down1)

    TAD_calling.subTAD_calling(tmp_dir1, up1, down1, aliase=aliases[0])

    tmp_dir2 = os.path.join(args.output, aliases[1])
    if not os.path.isdir(tmp_dir2):
        os.makedirs(tmp_dir2)

    up2, down2, files2 = TAD_calling.file_split(file[1], outfile=tmp_dir2, print_subcontact=0, aliases=aliases[1])

    print('auto:', up2, down2)

    if args.up_cutoff != 'None' and args.down_cutoff != 'None':
        up_cutoff = args.up_cutoff.split(',')
        up2 = float(up_cutoff[1])
        down_cutoff = args.down_cutoff.split(',')
        down2 = float(down_cutoff[1])
        print('manual:', up2, down2)

    TAD_calling.subTAD_calling(tmp_dir2, up2, down2, aliase=aliases[1])

    for i in range(len(files1)):
        contact_map_file1 = files1[i]
        contact_map_file2 = files2[i]
        TAD_file1 = contact_map_file1+'.TAD.txt'
        TAD_file2 = contact_map_file2 + '.TAD.txt'
        base1 = os.path.basename(contact_map_file1)
        base2 = os.path.basename(contact_map_file2)
        start1 = int(base1.split('.')[1])
        start2 = int(base2.split('.')[1])
        end = base2.split('.')[2]
        if start1 != start2:
            print('input file wrong')
            break

        TAD1 = np.loadtxt(TAD_file1)
        TAD2 = np.loadtxt(TAD_file2)

        try:
            TAD1[:, 1:3] = TAD1[:, 1:3]-start1
            TAD2[:, 1:3] = TAD2[:, 1:3]-start2
        except:
            pass

        map1 = np.loadtxt(contact_map_file1)
        map2 = np.loadtxt(contact_map_file2)

        try:
            TAD_split.write_split_tad(TAD1, TAD2, contact_map_file1, contact_map_file2, map1, map2, up1, up2,
                                      args.adjust_quality, args.output, start1, end, aliases[0], aliases[1])

            TAD_split.write_split_tad(TAD2, TAD1, contact_map_file2, contact_map_file1, map2, map1, up2, up1,
                                      args.adjust_quality, args.output, start1, end, aliases[1], aliases[0])

        except:
            pass

    a = TAD_split.merge_outputfile(aliases[0], aliases[1], args.output)
    if a[0].shape != 0:
        np.savetxt(os.path.join(args.output, aliases[0] + '->' + aliases[1] + '.all.merge.txt'), a[0])
    if a[1].shape != 0:
        np.savetxt(os.path.join(args.output, aliases[0] + '->' + aliases[1] + '.all.split.txt'), a[1])

    b = TAD_split.merge_outputfile(aliases[1], aliases[0], args.output)
    if b[0].shape != 0:
        np.savetxt(os.path.join(args.output, aliases[1] + '->' + aliases[0] + '.all.merge.txt'), b[0])
    if b[1].shape != 0:
        np.savetxt(os.path.join(args.output, aliases[1] + '->' + aliases[0] + '.all.split.txt'), b[1])




def TAD_similarity(command='TAD_similarity'):
    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython3 TADsplimer.py TAD_similarity <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python3 TADsplimer.py TAD_similarity -h\n")
        return 0

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython3 TADsplimer.py TAD_similarity <contact_map_file_paths> "
                                           "[optional arguments]\n\n", description='')
    parser.add_argument('command', default=None,
                        help="set as 'TAD_similarity' to calculate similarity of two TADs")

    parser.add_argument('-c', '--contact_maps', dest="contact_map", default=None,
                        help="paths to Hi-C contact maps in two conditions. paths must be separated by the comma ','.")

    parser.add_argument('-t', '--TAD', dest="TAD", default=None,
                        help="input files of TADs for two compared Hi-C contact maps. Paths must be separated by the"
                             " comma ','.")

    parser.add_argument('-o', '--output', dest="output", default=None,
                        help="path to output files")

    args = parser.parse_args()


    file1 = args.contact_map.split(',')
    file2 = args.TAD.split(',')
    map1 = np.loadtxt(file1[0])
    map2 = np.loadtxt(file1[1])
    TAD1 = np.loadtxt(file2[0])
    TAD2 = np.loadtxt(file2[1])
    head1, tail1 = os.path.split(file1[0])
    base1 = os.path.splitext(tail1)
    head2, tail2 = os.path.split(file1[1])
    base2 = os.path.splitext(tail2)

    TAD = np.concatenate((TAD1, TAD2), axis=0)
    scc = similarity_score.similarity_scc(map1, map2, TAD)
    Laplacian = similarity_score.similarity_Laplacian(map1, map2, TAD)
    hash = similarity_score.hash_similarity(map1, map2, TAD)
    np.savetxt(os.path.join(args.output, base1[0]+'_'+base2[0]+'_scc.txt'), scc, delimiter='\t')
    np.savetxt(os.path.join(args.output, base1[0]+'_'+base2[0]+'_laplacian.txt'), Laplacian, delimiter='\t')
    np.savetxt(os.path.join(args.output, base1[0]+'_'+base2[0]+'_hash.txt'), hash, delimiter='\t')


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == 'TAD_calculator':
            TAD_calculator(command='TAD_calculator')
        elif sys.argv[1] == 'split_TADs':
            split_TADs(command='split_TADs')
        elif sys.argv[1] == 'split_TADs_alternate':
            split_TADs_alternate(command='split_TADs_alternate')
        elif sys.argv[1] == 'TAD_similarity':
            TAD_similarity(command='TAD_similarity')
        else:
            printHelp()
    else:
        print('\nTADsplimer version 1.0.1')
        print('For a list of functions in TADsplimer, please try:\npython3 TADsplimer.py -h')
        print('')








