#!/share/app/python-3.4.3/bin/python3
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "spikeliu"

import sys
import argparse
import re

def parse_blat_file(blat_file):
    contig_in_ref = {}
    with open(blat_file, 'r') as f:
        # start and end pos of the contig on the ref, ignoring mismatch
        # end - start = length of the contig
        tstart = tend = 0
        for line in f:
            line = line.strip().split()
            if int(line[0]) < 0.95 * int(line[10]):
                continue
            if int(line[5]) + int(line[7]) > 20:
               #continue
               pass
            if line[8] == '+':
                tstart = int(line[15]) - int(line[11])
                tend = int(line[16]) + (int(line[10]) - int(line[12]))
            else:
                tstart = int(line[15]) - (int(line[10]) - int(line[12]))
                tend = int(line[16]) + int(line[11])
            # key = qname, value = (ref_name, strand, start, end)
            contig_in_ref[line[9]] = (line[13], line[8], tstart, tend)
    # print(contig_in_ref.__sizeof__())
    return contig_in_ref

def dynamic_search_contig(contig_in_scaff,
                          contig_in_ref,
                          min_overlap=0,
                          reverse=False):
    contig_cnt = len(contig_in_scaff)
    max_score = [-sys.maxsize - 1] * contig_cnt
    max_score_id = [-2] * contig_cnt
    for i in range(contig_cnt):
        contig_name_i, strand_i, start_i, length_i = contig_in_scaff[i]
        if contig_name_i not in contig_in_ref:
            msg = contig_name_i + '\t' + str(max_score[i])
            # print(msg)
            continue
        # detect inversion
        if reverse and strand_i == contig_in_ref[contig_name_i][1]:
            max_score_id[i] = -3
            continue
        elif not reverse and strand_i != contig_in_ref[contig_name_i][1]:
            max_score_id[i] = -3
            continue
        '''
        if contig_name_i in inversed_contig:
            max_score_id[i] = -3
            continue
        '''
        msg = contig_name_i + '\t' + str(length_i)
        max_score[i] = length_i
        max_score_id[i] = -1
        ref_name_i = contig_in_ref[contig_name_i][0]
        for j in range(i):
            contig_name_j, strand_j, start_j, length_j = contig_in_scaff[j]
            if contig_name_j not in contig_in_ref:
                continue
            if (ref_name_i != contig_in_ref[contig_name_j][0]):
                continue
            if reverse and strand_j == contig_in_ref[contig_name_j][1]:
                continue
            elif not reverse and strand_j != contig_in_ref[contig_name_j][1]:
                continue
            '''
            if contig_name_j in inversed_contig:
                continue
            '''
            contig_i = contig_in_ref[contig_name_i]
            contig_j = contig_in_ref[contig_name_j]
            end_i = start_i + length_i
            end_j = start_j + length_j

            # calculate 'fake' overlap
            shorter = max(length_i, length_j)
            bigger_start = max(contig_i[2], contig_j[2])
            smaller_end = min(contig_i[3], contig_j[3])
            # having overlap
            overlap = smaller_end - bigger_start
            duplicate = False
            # if overlap > min_overlap and overlap / float(shorter) >= 0.8:
            if (overlap - min_overlap) / float(shorter) >= 0.8:
                # print('overlap', overlap)
                duplicate = True

            score = 0
            # if not on the same ref, give a minimum int
            if contig_i[0] != contig_j[0]:
                score = - sys.maxsize - 1
            # if real duplicate
            elif duplicate:
                score = - sys.maxsize - 1
            else:
                # assume j starts behind i
                if reverse:
                    gap_in_ref = contig_j[2] - contig_i[3]
                # assume i starts behind j
                else:
                    gap_in_ref = contig_i[2] - contig_j[3]
                # print ('gap_in_ref', gap_in_ref)
                # always assume i starts behind j
                # if () or ():
                gap_in_scaff = start_i - (end_j)
                # print ('gap_in_scaff', gap_in_scaff)
                gap = abs(gap_in_scaff - gap_in_ref)
                # length = end_i - end_j
                # print ('gap_penalty', gap)

                if reverse:
                    length = -(contig_i[2] - contig_j[2])
                else:
                    length = contig_i[3] - contig_j[3]
                # print('length', length, i, j, reverse)
                if length < 0:
                    length = 0
                if length > length_i:
                    length = length_i
                # print('length', length, i, j)
                score =  max_score[j] + length - gap / 1.0
            msg += '\t' + str(score)
            if score >= max_score[i]:
                max_score[i] = score
                max_score_id[i] = j
        msg += '\t' + str(max_score[i])
        msg += '\t' + str(max_score_id[i])
        # print(msg)
    return max_score, max_score_id

def detect_inverson(contig_in_scaff, contig_in_ref):
    same_strand = set()
    diff_strand = set()
    same_length = diff_length = 0
    for i in range(len(contig_in_scaff)):
        contig_name, strand, start, length = contig_in_scaff[i]
        if contig_name not in contig_in_ref:
            continue
        strand_on_ref = contig_in_ref[contig_name][1]
        if strand == strand_on_ref:
            same_strand.add(contig_name)
            same_length += length
        else:
            diff_strand.add(contig_name)
            diff_length += length

    # suppose same_stand has most members
    # diff_strand holds the contig inversed
    if len(same_strand) < len(diff_strand):
        same_strand, diff_strand = diff_strand, same_strand
    elif len(same_strand) == len(diff_strand):
        if same_length < diff_length:
            same_strand, diff_strand = diff_strand, same_strand
    return diff_strand

def evaluate_scaff_struct(scaff_name, scaff_length, contig_in_ref, contig_in_scaff,
                          min_overlap):
    # print('>' + scaff_name)
    msg = '>' + scaff_name  + '\n'

    # detect inversion
    # inversed_contig = detect_inverson(contig_in_scaff, contig_in_ref)
    # inversed_contig = set()

    # print('[FORWARD]')
    max_score, max_score_id = dynamic_search_contig(contig_in_scaff,
                                                    contig_in_ref,
                                                    min_overlap,
                                                    False)

    # print('[BACKWARD]')
    max_score_reverse, max_score_id_reverse = dynamic_search_contig(
            contig_in_scaff, contig_in_ref, min_overlap, True)

    if max(max_score_reverse) > max(max_score):
        max_score = max_score_reverse
        max_score_id = max_score_id_reverse
    # trace back
    '''
    if sum(max_score_id) == -len(max_score_id) * 2:
        print(scaff_name + '\'s contigs are all not on the ref')
        return False
    '''
    contig_pass = [False] * len(max_score)
    max_id = 0
    score = max_score[0]
    for i in range(1, len(max_score)):
        if max_score[i] >= score:
            max_id = i
            score = max_score[i]
    if max_score_id[max_id] != -2 and max_score_id[max_id] != -3:
        # first contig has been found
        # contig_pass.append(contig_in_scaff[max_id][0])
        contig_pass[max_id] = True
        contig_id = max_score_id[max_id]
        while(1):
            if contig_id == -1:
                break
            # contig_pass.append(contig_in_scaff[contig_id][0])
            contig_pass[contig_id] = True
            contig_id = max_score_id[contig_id]

    # print contigs with different status
    pass_cnt = mis_cnt = inver_cnt = trans_cnt = 0
    mismatch = inversion = translocation = 0
    ref_winner = ""
    contig_name_lead = contig_in_scaff[max_id][0]
    if contig_name_lead in contig_in_ref:
        ref_winner = contig_in_ref[contig_name_lead][0]

    for i in range(len(contig_in_scaff)):
        contig_name, strand, start, length = contig_in_scaff[i]
        msg += contig_name
        if  contig_pass[i]:
            pass_cnt += 1
            msg += ' (chosen)'
        elif max_score_id[i] == -2:
            msg += ' (not on ref)'
        elif contig_in_ref[contig_name][0] != ref_winner:
            msg += ' (translocated)'
            trans_cnt += 1
            translocation += length
        elif max_score_id[i] == -3:
            msg += ' (inversed)'
            inver_cnt += 1
            inversion += length
        else:
            msg += ' (ignored)'
            # caution
            mis_cnt += 1
            mismatch += length
        msg += ' ' + str(start) + ' ' + strand + ' ' + str(length)
        if contig_name in contig_in_ref:
            ref_name, strand, start, end = contig_in_ref[contig_name]
            msg += ' info on ref: ' + ref_name + ' ' + strand + ' ' \
                       + str(start) + ' ' + str(end)
        else:
            msg += ' not available'
        msg += '\n'
        # print(msg)

    error_length = mismatch + inversion + translocation
    iden = 1 - error_length / float(scaff_length)
    # print mismatch length
    '''
    print('# of chosen contig:', pass_cnt, ', ratio:',
          float(pass_cnt) / len(contig_in_scaff), '\nidentity:',
          iden)
    '''
    ratio = float(pass_cnt) / len(contig_in_scaff)
    msg += '[STAT] # of chosen contig: ' + str(pass_cnt) + ', ratio: ' + str(ratio) + \
            '\n[STAT] identity: ' + str(iden) + '\n'
    '''
    if scaff_length >= 1000000 and iden < 0.99:
        print('[DEBUG] %d %f' % (scaff_length, iden))
    '''

    msg += '[STAT] error length: ' + str(error_length)
    if (error_length > 400000):
        msg += ' (big mistake)'
    msg += '\n'

    if iden >= 0.95:
        # print('[SUCCESS]')
        msg += '[SUCCESS]'
        print(msg)
        return True, mis_cnt, mismatch, inver_cnt, inversion, trans_cnt, translocation
    else:
        # print('[FAIL]')
        msg += '[FAIL]'
        print(msg)
        return False, mis_cnt, mismatch, inver_cnt, inversion, trans_cnt, translocation

def parse_scaff_file(scaff_file, contig_in_ref, min_overlap, get_length_from_name=True):
    success = cnt = 0
    mismatch = inversion = translocation = 0
    total_length = 0
    mismatch_cnt = inversion_cnt = translocation_cnt = 0
    with open(scaff_file, 'r') as f:
        contig_in_scaff = []
        line = f.readline().strip()
        scaff_name = line[1:]
        scaff_length = 0
        if get_length_from_name:
            scaff_length = int(line.split()[-1])
        for line in f:
            line = line.strip()
            if (line[0] == '>'):
                cnt += 1
                if len(contig_in_scaff) > 0:
                    if not get_length_from_name:
                        scaff_length = contig_in_scaff[-1][2] + \
                                       contig_in_scaff[-1][3]

                    is_succ, mis_cnt, mis, inver_cnt, inver, t_cnt, t = \
                                evaluate_scaff_struct(scaff_name,
                                                      scaff_length,
                                                      contig_in_ref,
                                                      contig_in_scaff,
                                                      min_overlap)
                    success += is_succ
                    mismatch_cnt += mis_cnt
                    mismatch += mis
                    inversion_cnt += inver_cnt
                    inversion += inver
                    translocation_cnt += t_cnt
                    translocation += t
                    total_length += scaff_length
                scaff_name = line[1:]
                if get_length_from_name:
                    scaff_length = int(line.split()[-1])
                contig_in_scaff = []
            else:
                line = line.split()
                # qname, strand, start, length
                contig_in_scaff.append((line[0], line[2],
                                        int(line[1]), int(line[3])))
        if len(contig_in_scaff) > 0:
            if not get_length_from_name:
                scaff_length = contig_in_scaff[-1][2] + \
                        contig_in_scaff[-1][3]

            is_succ, mis_cnt, mis, inver_cnt, inver, t_cnt, t = \
                                   evaluate_scaff_struct(scaff_name,
                                                         scaff_length,
                                                         contig_in_ref,
                                                         contig_in_scaff,
                                                         min_overlap)
            success += is_succ
            mismatch_cnt += mis_cnt
            mismatch += mis
            inversion_cnt += inver_cnt
            inversion += inver
            translocation_cnt += t_cnt
            translocation += t
            total_length += scaff_length
    print('[SUMMARY] # of passed scaff:', success, 'total:', cnt,
          'ratio:', float(success) / cnt)
    inconsistency=translocation+inversion+mismatch
    print('[SUMMARY] # of insertion or relocation: %d' % mismatch_cnt,
          '\n[SUMMARY] length of insertion or relocation: %d bp' % mismatch,
          '\n[SUMMARY] # of inversion: %d' % inversion_cnt,
          '\n[SUMMARY] length of inversion: %d bp' % inversion,
          '\n[SUMMARY] # of translocation: %d' % translocation_cnt,
	  '\n[SUMMARY] length of translocation: %d bp' % translocation,
          '\n[SUMMARY] length of inconsistencies: %d bp' % inconsistency,
          '\n[SUMMARY] total length: %d bp' % total_length)

# not used
def evaluate_scaff(scaff_name, contig_in_scaff):
    for i in range(len(contig_in_scaff)-1):
        contig_name_i, strand_i, start_i, length_i = contig_in_scaff[i]
        contig_name_j, strand_j, start_j, length_j = contig_in_scaff[i+1]
        if start_i >= start_j or start_i + length_i >= start_j + length_j:
            print(scaff_name, contig_name_i, contig_name_j)

# not used
def check_gap(scaff_file):
    fail_cnt = 0
    with open(scaff_file, 'r') as f:
        contig_in_scaff = []
        scaff_name = f.readline().strip().split()[0][1:]
        for line in f:
            if (line[0] == '>'):
                evaluate_scaff(scaff_name,
                               contig_in_scaff)
                scaff_name = line.strip().split()[0][1:]
                contig_in_scaff = []
            else:
                line = line.strip().split()
                # qname, strand, start, length
                contig_in_scaff.append((line[0], line[2],
                                        int(line[1]), int(line[3])))
        evaluate_scaff(scaff_name,
                       contig_in_scaff)

# check gap difference based on evaluation result.
def gap_diff(eval_file, out_file):
    with open(eval_file, 'r') as f:
        # [(start_r, end_r, start_s, end_s), ...]
        gap_differences = []
        gaps = []
        for line in f:
            line = line.strip()
            #print(line)
            if line.startswith('>'):
                if len(gaps) >= 2:
                    for i in range(len(gaps) - 1):
                        gap_r = gaps[i + 1][0] - gaps[i][1]
                        gap_s = gaps[i + 1][2] - gaps[i][3]
                        gap_differences.append(gap_r - gap_s)
                gaps = []
                continue
            if re.match(r'\[|#|identity', line):
                continue
            content = line.split()
            if content[1].startswith('(ignored)') or content[1].startswith('(inversed)'):
                continue
            start_s = int(content[2])
            end_s = start_s + int(content[4])
            start_r = int(content[10])
            end_r = int(content[11])
            gap = (start_r, end_r, start_s, end_s)
            gaps.append(gap)
    # deal with last scaffold
    if len(gaps) >= 2:
        for i in range(len(gaps) - 1):
            gap_r = abs(gaps[i + 1][0] - gaps[i][1])
            gap_s = abs(gaps[i + 1][2] - gaps[i][3])
            gap_differences.append(gap_r - gap_s)
    with open(out_file, 'w') as f:
        for diff in gap_differences:
            f.write(str(diff) + '\n')

def main():
    # using argparse to parse the command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("scaff_file", help="scaffold input", type = str)
    parser.add_argument("blat_file", help="blat input", type = str)
    parser.add_argument("--min-overlap",
                        help="minimum overlap, could be Kmer size",
                        metavar='min_overlap', default = 0, type = int)
    parser.add_argument("--get-length-from-name",
                        help="get scaff length from its name, defaull is True",
                        metavar='min_overlap', default = True, type = bool)
    args = parser.parse_args()
    contig_in_ref = parse_blat_file(args.blat_file)
    parse_scaff_file(args.scaff_file, contig_in_ref, args.min_overlap, args.get_length_from_name)

if __name__ == '__main__':
    main()

