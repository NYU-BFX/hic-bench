import numpy as np
from itertools import combinations
import similarity_score
import os


def TAD_matching(TAD1, TAD2):
    """
    matching TADs in two conditions
    :param TAD1: TAD in the first condition
    :param TAD2: TAD in the second condition
    :return: TADs in two conditions, TADs only in the first condition, TADs only in the second condition
    """
    flag1 = []
    flag2 = []
    matched_TAD = np.array([0, 0, 0])
    try:
        for i in range(0, TAD1.shape[0]):
            for j in range(0, TAD2.shape[0]):
                if abs(TAD1[i, 2] - TAD2[j, 2]) <= 10 and abs(TAD1[i, 1] - TAD2[j, 1]) <= 10:
                    matched_TAD = np.vstack([matched_TAD, (TAD1[i, :] + TAD2[j, :]) / 2])
                    flag1.append(i)
                    flag2.append(j)
        TAD1_only = np.delete(TAD1, flag1, axis=0)
        TAD2_only = np.delete(TAD2, flag2, axis=0)
    except:
        try:
            TAD2 = TAD2.reshape([1, 3])
            for i in range(0, TAD1.shape[0]):
                for j in range(0, TAD2.shape[0]):
                    if abs(TAD1[i, 2] - TAD2[j, 2]) <= 10 and abs(TAD1[i, 1] - TAD2[j, 1]) <= 10:
                        matched_TAD = np.vstack([matched_TAD, (TAD1[i, :] + TAD2[j, :]) / 2])
                        flag1.append(i)
                        flag2.append(j)
            TAD1_only = np.delete(TAD1, flag1, axis=0)
            TAD2_only = np.delete(TAD2, flag2, axis=0)
        except:
            TAD1 = TAD1.reshape([1, 3])
            for i in range(0, TAD1.shape[0]):
                for j in range(0, TAD2.shape[0]):
                    if abs(TAD1[i, 2] - TAD2[j, 2]) <= 10 and abs(TAD1[i, 1] - TAD2[j, 1]) <= 10:
                        matched_TAD = np.vstack([matched_TAD, (TAD1[i, :] + TAD2[j, :]) / 2])
                        flag1.append(i)
                        flag2.append(j)
            TAD1_only = np.delete(TAD1, flag1, axis=0)
            TAD2_only = np.delete(TAD2, flag2, axis=0)
    try:
        return matched_TAD[1:, :], TAD1_only, TAD2_only
    except:
        return matched_TAD, TAD1_only, TAD2_only


def get_matched_TAD_dictionary(range1, range2, TAD_s):

    def RangInList(L, TAD):
        flag1 = 1
        try:
            for t in TAD:
                if abs(L[1] - t[1]) < 10 and abs(L[2] - t[2]) < 10:
                    flag1 = 0
        except:
            flag1 = 1
        return flag1

    flag = 0
    for t in range2:
        if abs(t[1] - range1[1]) < 10:
            flag = flag + 1
        if abs(t[2] - range1[2]) < 10:
            flag = flag + 1
    if flag >= 2:
        matched_TAD_dictionary = []
        for i in range(2, range2.shape[0] + 1):
            TAD_combination = list(combinations(range(range2.shape[0]), i))
            for j, t in enumerate(TAD_combination):
                x = range2[t, :]
                length = 0
                flag3 = []
                for s in x:
                    length = length + s[2] - s[1]
                    f2 = RangInList(s, TAD_s.tolist())
                    flag3.append(f2)
                r = length / (range1[2] - range1[1])
                if 0.8 < r < 1.1:
                    flag2 = 0
                    for l in range(x.shape[0] - 1):
                        if x[l + 1, 1] - x[l, 2] > 20 or x[l, 2] - x[l + 1, 1] > 20:
                            flag2 = 1
                    if flag2 == 0:
                        matched_TAD_dictionary.append([x, r])
        return matched_TAD_dictionary


def corner_split_score(TAD1, TAD2, TAD_s):
    '''
    calculate the corner split score for two set of TADs
    :param TAD1: all TADs for calculation
    :param TAD2: divide region
    :param TAD_s: similar region
    :return: merged_TADs: merged TADs, D2: split TADs with coverage rate, split_TAD: split TADs without coverage rate
    '''
    
    matched_TAD = np.empty(shape=[0, 3])
    merged_TAD = []
    split_TAD_and_CSR = []
    split_TAD = []
    count = 0
    try:
        for i, t1 in enumerate(TAD1):
            for t2 in TAD2:
                if t1[1] - t2[1] < 20 and t2[2] - t1[2] < 20:
                    matched_TAD = np.append(matched_TAD, [t2], axis=0)
            D = get_matched_TAD_dictionary(t1, matched_TAD, TAD_s)

            if D is not None and len(D) != 0:
                m = []
                for j in range(0, len(D)):
                    m.append(D[j][0].shape[0])
                for j in range(0, len(D)):
                    merged_TAD.append(t1)
                    split_TAD_and_CSR.append(D[j])
                    split_TAD.append([D[j][0][:, 1:3], count])
                    count = count + 1
            matched_TAD = np.empty(shape=[0, 3])
    except:
        try:
            TAD1 = TAD1.reshape([1, 3])
            for i, t1 in enumerate(TAD1):
                for t2 in TAD2:
                    if t1[1] - t2[1] < 20 and t2[2] - t1[2] < 20:
                        matched_TAD = np.append(matched_TAD, [t2], axis=0)
                D = get_matched_TAD_dictionary(t1, matched_TAD, TAD_s)
                if D is not None and len(D) != 0:
                    m = []
                    for j in range(0, len(D)):
                        m.append(D[j][0].shape[0])
                    for j in range(0, len(D)):
                        if D[j][0].shape[0] == max(m):
                            merged_TAD.append(t1)
                            split_TAD_and_CSR.append(D[j])
                            split_TAD.append([D[j][0][:, 1:3], count])
                            count = count + 1
                matched_TAD = np.empty(shape=[0, 3])
        except:
            TAD2 = TAD2.reshape([1, 3])
            for i, t1 in enumerate(TAD1):
                for t2 in TAD2:
                    if t1[1] - t2[1] < 20 and t2[2] - t1[2] < 20:
                        matched_TAD = np.append(matched_TAD, [t2], axis=0)
                D = get_matched_TAD_dictionary(t1, matched_TAD, TAD_s)
                if D is not None and len(D) != 0:
                    m = []
                    for j in range(0, len(D)):
                        m.append(D[j][0].shape[0])

                    for j in range(0, len(D)):
                        if D[j][0].shape[0] == max(m):
                            merged_TAD.append(t1)
                            split_TAD_and_CSR.append(D[j])
                            split_TAD.append([D[j][0][:, 1:3], count])
                            count = count + 1
                matched_TAD = np.empty(shape=[0, 3])
    return np.array(merged_TAD), np.array(split_TAD_and_CSR), np.array(split_TAD)


def average_interaction(split_TAD, map1, map2, p1, p2):
    '''
    :param split_TAD: split TADs
    :param map1: contact map in the first condition
    :param map2: contact map in the second condition
    :return: average interaction between two split TADs
    '''
    average_interaction1 = []
    average_interaction2 = []
    average_diff = []
    for count, d0 in enumerate(split_TAD):
        d = d0[0]
        avr_all1_temp = np.empty(shape=[d.shape[0], d.shape[0]])
        avr_all2_temp = np.empty(shape=[d.shape[0], d.shape[0]])
        for i in range(d.shape[0]):
            for j in range(d.shape[0]):
                m1 = map1[int(d[i][0]): int(d[i][1]), int(d[j][0]): int(d[j][1])]
                m2 = map2[int(d[i][0]): int(d[i][1]), int(d[j][0]): int(d[j][1])]
                avr1 = np.mean(m1[m1 < p1])
                avr_all1_temp[i, j] = avr1
                avr2 = np.mean(m2[m2 < p2])
                avr_all2_temp[i, j] = avr2
        average_interaction1.append([avr_all1_temp, count])
        average_interaction2.append([avr_all2_temp, count])
        average_diff.append([avr_all1_temp - avr_all2_temp, count])
    return average_interaction1, average_interaction2, average_diff


def split_region(f1, f2, split_TAD, merged_TADs, p1, p2, fold):
    '''
    main function for divided region detection
    :param f1: file for contact map in the first condition (split)
    :param f2: file for contact map in the second condition
    :param split_TAD: divided region in cell line 1
    :param merged_TADs: merged TADs in the second condition
    :param TAD_s: matched TADs
    :return: divid1: mean interaction in divide region, divid2: mean interaction in union region, split_TAD_location: location of
    divided region, merged_TAD_location: location of union region, dlt: similar region after removing differential region
    '''

    def get_interaction(idx, avr):
        '''
        calculate average interaction
        :param idx: index of real region
        :param avr: interact of all region
        :return: average interaction
        '''
        inter = np.empty(shape=(0, 3))
        for t in avr:
            if t[1] in idx:
                T = t[0]
                for i in range(T.shape[0] - 1):
                    inter_tmp = [T[i, i], T[i, i + 1], T[i + 1, i + 1]]
                    inter = np.append(inter, [inter_tmp], axis=0)
        return inter

    def get_corner_ratio(avr, fold):
        '''
        the ratio of interaction in conner region/ interaction in diagonal region
        :param avr: mean of interaction in each region
        :return: matrix of ratio, flag == 1 is real merge region
        '''

        def corner_split_ratio(contact_map):
            R = np.zeros((contact_map.shape[0] - 1, contact_map.shape[0] - 1))
            for i in range(contact_map.shape[0] - 1):
                for j in range(contact_map.shape[0] - 1):
                    R[i, j] = (contact_map[i, j + 1] * 2) / (contact_map[i, j] + contact_map[i + 1, j + 1])
            return R

        ratio = []
        for count, t in enumerate(avr):
            contact_map = t[0]
            R = corner_split_ratio(contact_map)
            x1 = R.diagonal()[R.diagonal() > 0.47].shape[0]
            x2 = R.diagonal()[R.diagonal() < 0.47].shape[0]
            flag = 0
            if R.shape[0] == 1:
                if R[0, 0] > 0.45:
                    flag = 1
            else:
                if x2 < x1:
                    flag = 1
            ratio.append([R / fold, count, flag])
        return ratio

    def get_split_region(ratio1, ratio2):
        '''
        get split site
        :param ratio1: split region
        :param ratio2: merge region
        :return: real split site
        '''
        idx = []
        ratio_diff = []
        for i in range(len(ratio1)):
            if ratio1[i][2] != 1:
                ratio_tmp = ratio2[i][0] - ratio1[i][0]
                if np.mean(ratio_tmp.diagonal()) > 0.09:  # real divide region by cutoff: -0.1
                    idx.append(ratio2[i][1])
                    ratio_diff.append(np.mean(ratio_tmp.diagonal()))
        return idx, ratio_diff

    def get_split_TAD(split_TAD, idx):
        '''
        get the split TADs
        :param split_TAD: all divided region
        :param idx: real region
        :return: location of real region
        '''
        reg = np.empty(shape=(0, 2))
        for d in split_TAD:
            if d[1] in idx:
                reg = np.append(reg, d[0], axis=0)
        return reg
    
    map1 = np.loadtxt(f1)
    map2 = np.loadtxt(f2)
    avr1, avr2, average_diff = average_interaction(split_TAD, map1, map2, p1,
                                                 p2)
    ratio1 = get_corner_ratio(avr1, 1)
    ratio2 = get_corner_ratio(avr2, fold)
    idx, ratio_diff = get_split_region(ratio1, ratio2)
    ratio_diff = np.array(ratio_diff)

    if len(idx) > 0:
        divid1 = get_interaction(idx, avr1)
        divid2 = get_interaction(idx, avr2)
        merged_TAD_location = merged_TADs[idx, :]
        merged_TAD_location = np.c_[merged_TAD_location, ratio_diff]
        split_TAD_location = get_split_TAD(split_TAD, idx)
        split_TAD_location = np.insert(split_TAD_location, 0, np.array(range(split_TAD_location.shape[0])).transpose(),
                                       axis=1)
        return divid1, divid2, split_TAD_location, merged_TAD_location
    else:
        return 0, 0, 0, 0


def find_t1(TAD):
    def is_t1(t, TAD):
        '''
        t is T1 TAD return 1, else return 0
        '''
        flag = 1
        for t0 in TAD:
            if not np.array_equal(t0, t):
                if t[1] - t0[1] <= 0 and t0[2] - t[2] <= 0:
                    flag = 0
        return flag

    t1 = np.empty(shape=[0, 3])
    for t in TAD:
        if is_t1(t, TAD) == 1:
            t1 = np.append(t1, [t], axis=0)
    return t1


def get_ratio(map, TAD, up):
    def neighbor_TAD(TAD):
        '''
        find each pair of nearby TAD
        '''
        pair = []
        for i in range(TAD.shape[0]):
            for j in range(TAD.shape[0]):
                if abs(TAD[i, 2] - TAD[j, 1]) < 20:
                    pair.append([TAD[i], TAD[j]])
        return np.array(pair)

    def calculate_ratio(map, pair):
        ratio = []
        for T in pair:
            t1 = T[0]
            t2 = T[1]
            m1 = map[int(t1[1]):int(t1[2]), int(t1[1]):int(t1[2])]
            m2 = map[int(t2[1]):int(t2[2]), int(t2[1]):int(t2[2])]
            m3 = map[int(t1[1]):int(t1[2]), int(t2[1]):int(t2[2])]
            mean_m1 = np.mean(m1[m1 < up])
            mean_m2 = np.mean(m2[m2 < up])
            mean_m3 = np.mean(m3[m3 < up])
            r = mean_m3 / (mean_m1 + mean_m2)
            ratio.append(r)
        return ratio

    pair = neighbor_TAD(TAD)
    ratio = calculate_ratio(map, pair)
    return ratio


def write_split_tad(TAD1, TAD2, contact_map_file1, contact_map_file2, map1, map2, up1, up2,
                    adjust_quality, output, start1, end, aliases1, aliases2):

    TAD_s, TAD1_only, TAD2_only = TAD_matching(TAD2, TAD1)
    merged_TADs, D2, split_TAD = corner_split_score(TAD1, TAD2, TAD_s)

    if adjust_quality == 0:
        split1_1, split2_1, split_TAD_location, merged_TAD_location = split_region(contact_map_file2,
                                                                contact_map_file1, split_TAD, merged_TADs, up2, up1, 1)
    elif adjust_quality == 1:
        ratio1 = get_ratio(map1, TAD1, up1)
        ratio2 = get_ratio(map2, TAD2, up2)
        fold = np.mean(ratio1) / np.mean(ratio2)
        split1_1, split2_1, split_TAD_location, merged_TAD_location = split_region(contact_map_file2, contact_map_file1,
                                                                  split_TAD, merged_TADs, up2, up1, fold)
    else:
        print("\nusage:\npython3 DACTAD.py TAD_calculator <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python DACTAD.py TAD_calculator -h\n")
    try:
        if split1_1 == 0:
            print('No split TAD')
    except:
        if (not os.path.exists(output)):
            print('path does not exist\n')
        else:
            scc = similarity_score.similarity_scc(map1, map2, merged_TAD_location[:, 0:3])
            Laplacian = similarity_score.similarity_Laplacian(map1, map2, merged_TAD_location[:, 0:3])
            hash = similarity_score.hash_similarity(map1, map2, merged_TAD_location[:, 0:3])
            loc_u2 = np.c_[merged_TAD_location, scc[:, 3], Laplacian[:, 3], hash[:, 3]]
            np.savetxt(os.path.join(output, aliases1 + '->' + aliases2 + '.' + str(start1) + '.' + end
                                    + '.split.txt'), split_TAD_location)
            np.savetxt(os.path.join(output, aliases1 + '->' + aliases2 + '.' + str(start1) + '.' + end
                                    + '.merge.txt'), loc_u2)


def merge_outputfile(aliases1, aliases2, path):

    def TAD_redundance(t, TAD):
        '''
        flag == 0 to remove; flag == 1 not to remove
        '''
        flag = 1
        try:
            for tad in TAD:
                if abs(t[1] - tad[1]) < 20 and abs(t[2] - tad[2]) < 20:
                    flag = 0
        except: pass
        return flag

    def remove_redundance(TAD):
        idx = []
        for i in range(TAD.shape[0]):
            remove = TAD_redundance(TAD[i, :], TAD[(i+1):, :])
            if remove == 1:
                idx.append(i)
        TAD_nonredundance = (TAD[idx, :])
        TAD_nonredundance[:, 0] = range(TAD_nonredundance.shape[0])
        return TAD_nonredundance

    files = os.listdir(path)
    split_region = np.empty([1, 3])
    merge_region = np.empty([1, 7])
    for f in files:
        if aliases1+'->'+aliases2 in f:
            start_loc = f.split('.')[-4]
            TAD = np.loadtxt(os.path.join(path, f))
            if f.split('.')[-2] == 'split':
                try:
                    TAD = np.reshape(TAD, (1, 3))
                except:
                    pass
                TAD[:, 1:3] = TAD[:, 1:3] + int(start_loc)
                split_region = np.vstack([split_region, TAD])
            else:
                try:
                    TAD = np.reshape(TAD, (1, 7))
                except:
                    pass
                TAD[:, 1:3] = TAD[:, 1:3] + int(start_loc)
                merge_region = np.vstack([merge_region, TAD])

    split_region = split_region[1:, :]
    merge_region = merge_region[1:, :]
    split_region = split_region[np.argsort(split_region[:, 1])]
    merge_region = merge_region[np.argsort(merge_region[:, 1])]
    return (remove_redundance(merge_region), remove_redundance(split_region))
