#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "spikeliu"

import argparse
import gzip
import pysam
import re
import os
import shutil
import subprocess

class corrdinate:
    def __init__(self, qstart, qend, tstart, tend):
        self.qstart = qstart
        self.qend = qend
        self.tstart = tstart
        self.tend = tend

    def __eq__(self, other):
        return ((self.tstart, self.tend) == (other.tstart, other.tend))

    def __ne__(self, other):
        return ((self.tstart, self.tend) != (other.tstart, other.tend))

    def __lt__(self, other):
        return ((self.tstart, self.tend) < (other.tstart, other.tend))

    def __le__(self, other):
        return ((self.tstart, self.tend) <= (other.tstart, other.tend))

    def __gt__(self, other):
        return ((self.tstart, self.tend) > (other.tstart, other.tend))

    def __ge__(self, other):
        return ((self.tstart, self.tend) >= (other.tstart, other.tend))

    def set_all(self, qstart, qend, tstart, tend):
        self.qstart = qstart
        self.qend = qend
        self.tstart = tstart
        self.tend = tend

    def get_all(self):
        return self.qstart, self.qend, self.tstart, self.tend

def create_seq_map(fasta_in, stat=False):
    seq_map = {}
    length_list = []
    with pysam.FastxFile(fasta_in) as fin:
        for entry in fin:
            seq_map[entry.name] = entry.sequence.upper()
            if stat:
                length_list.append((len(entry.sequence), entry.name))

    if stat:
        print('seq stats')
        length_list.sort()
        for slen, sname in length_list:
            print(sname, slen)

    return seq_map



def rc(seq):
    bases = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'a':'t', 't':'a',
             'c':'g', 'g':'c', 'n':'n'}
    seq_re = seq[::-1]
    seq_rc = ''
    for b in seq_re:
        seq_rc += bases[b]
    return seq_rc

def create_coordinate(record):
    qstart = int(record[2])
    qend = qstart + int(record[4])
    tstart = int(record[10])
    tend = int(record[11])
    return corrdinate(qstart, qend, tstart, tend);

def deal_chosen(record, qseq, tseq):
    corrd = create_coordinate(record)
    qstart, qend, tstart, tend = corrd.get_all()
    '''
    # check whether 1st and last base on query are same to target
    if record[3] == '+' and qseq[qstart] != tseq[tstart]:
        print('wtf:', record[8], qstart, qseq[qstart], tstart, tseq[tstart])
        print(qseq[qstart:qend])
        print(tseq[tstart:tend])
    if record[3] == '+' and qseq[qend - 1] != tseq[tend - 1]:
        print('wtf:', record[8], qend, qseq[qend - 1], tend, tseq[tend - 1])
        print(qseq[qstart:qend])
        print(tseq[tstart:tend])
    '''
    '''
    while qstart < qend and tstart < tend:
        if qseq[qstart] == tseq[tstart]:
            break
        qstart += 1
        tstart += 1
    while qend > qstart and tend > tstart:
        if qseq[qend - 1] == tseq[tend - 1]:
            break
        qend -= 1
        tend -= 1
    '''
    if qstart == qend or tstart == tend:
        return None
    if record[3] == '-':
        return [corrdinate(-1, -1, tstart, tend), rc(qseq[qstart:qend])]
    return [corrdinate(-1, -1, tstart, tend), qseq[qstart:qend]]

def deal_inversed(record, qseq, tseq):
    return deal_chosen(record, qseq, tseq)
    #return None
    

def deal_translocated(record, qseq, tseq):
    return deal_chosen(record, qseq, tseq)
    return None

def deal_ignored(record, qseq, tseq):
    #return deal_chosen(record, qseq, tseq)
    return None

def deal_not_on_ref(record, qseq, tseq):
    return None

def parse_record(record, qseq, tseq):
    signal = record[1].strip('()')
    if signal == 'chosen':
        return deal_chosen(record, qseq, tseq)
    elif signal == 'inversed':
        return deal_inversed(record, qseq, tseq)
    elif signal == 'translocated':
        return deal_translocated(record, qseq, tseq)
    elif signal == 'ignored':
        return deal_ignored(record, qseq, tseq)
    elif signal == 'not_on_ref':
        return deal_not_on_ref(record, qseq, tseq)
    else:
        print('unkonw signal')
        exit(1)

def parse_eval_file(eval_in, query_seq_map, target_seq_map):
    replacements = {}
    with open(eval_in, 'r') as fin:
        line = fin.readline()
        if not line.startswith('>'):
            print('something is wrong with evaluation file, please check')
            exit(1)
        qname = line.strip().split()[0][1:]
        qseq = query_seq_map[qname]
        prev = None
        for line in fin:
            if line.startswith('>'):
                qname = line.strip().split()[0][1:]
                qseq = query_seq_map[qname]
                prev = None
                continue
            if line.startswith('['):
                continue
            record = line.strip().split()
            tname = record[8]
            tseq = target_seq_map[tname]
            res = parse_record(record, qseq, tseq)
            if not res:
                continue
            if tname in replacements:
                replacements[tname].append(res)
            else:
                replacements[tname] = [res]
            # extra check for ovalap
            if prev:
                pass
                #_, _, ts1, te1 = prev[0].get_all()
                #_, _, ts2, te2 = res[0].get_all()
            prev = res
                

    for _, replacement in replacements.items():
        replacement.sort()

    return replacements

def parse_cs_tag(qseq, cs_tag, max_length_to_replace):
    cs_tag = cs_tag[5:]
    p = re.compile(r'[:*+\-~]')
    signal = p.findall(cs_tag)
    content = p.split(cs_tag)[1:]
    if len(signal) != len(content):
        print('something wrong with cs tag')
        exit(1)
    i = 0
    qseq_new = ''
    off = 0
    while i < len(signal):
        if signal[i] == ':':
            end_pos = int(content[i])
            qseq_new += qseq[off:off + end_pos]
            off += end_pos
        elif signal[i] == '*':
            if content[i][1] == 'n':
                qseq_new += content[i][0].upper()
            else:
                qseq_new += content[i][1].upper()
            off += 1
        elif signal[i] == '+':
            if not max_length_to_replace or len(content[i]) <= max_length_to_replace:
                if content[i].count('n'):
                    qseq_new += content[i].replace('n', '').upper()
                else:
                    qseq_new += content[i].upper()
            end_pos = len(content[i])
            off += end_pos
        elif signal[i] == '-':
            if max_length_to_replace and len(content[i]) > max_length_to_replace:
                qseq_new += content[i].upper()
        elif signal[i] == '~':
            print('there should not be a tilde signal')
            exit(1)
        i += 1
    if off != len(qseq):
        print('weird, they should be the same')
        exit(1)
    if qseq_new.count('N') or qseq_new.count('n'):
        print('replace N failed')
        exit(1)
    return qseq_new

def rearrange_replacements(replacements, target_seq_map, tmp_path,
                           max_length_to_replace, num_thread):
    command = ' /home/wupei1/wupei1/software/minimap2 -c -x asm5 ' \
              ' --mask-level 0.9' \
              ' --min-occ 200 -g 2500 --score-N 2 --cs -t 1 '
    repl_length = 0
    parallel_command = ' /home/liuchao3/data/software/parallel-20200322/src/parallel --lb '
    fail_cnt = 0
    if not os.path.exists(tmp_path):
        if not os.path.isdir(tmp_path):
            os.makedirs(tmp_path)
    for tname, replacement in replacements.items():
        #if tname != 'ctg210':
        #    continue
        print('start to process', tname, len(replacement))
        i = 0
        while i < len(replacement) - 1:
            corrd1, _ = replacement[i]
            corrd2, _ = replacement[i + 1]
            _, _, tstart1, tend1 = corrd1.get_all()
            _, _, tstart2, tend2 = corrd2.get_all()
            # overlapping
            if tstart2 < tend1:
                # if one contains another
                if tend2 <= tend1:
                    replacement.pop(i + 1)
                    continue
                elif tstart1 == tstart2 and tend1 <= tend2:
                    replacement.pop(i)
                    continue
                else: # simple overlap
                    replacement[i][0].set_all(-1, -1, tstart1, tstart2)
                    i += 1
            else:
                i += 1
        print('after get rid of overlap', tname, len(replacement))
        tmp_path_para = tmp_path + '/' + tname
        if not os.path.exists(tmp_path_para):
            if not os.path.isdir(tmp_path_para):
                os.makedirs(tmp_path_para)
        
        process_num = 0
        # remap qseq to tseq                
        for i in range(len(replacement)):
            corrd, qseq = replacement[i]
            _, _, tstart, tend = corrd.get_all() 
            tseq = target_seq_map[tname][tstart:tend]
            with open(tmp_path_para + '/' + tname + '_q_' + str(i) + '.fa', 'w') as qf, \
                 open(tmp_path_para + '/' + tname + '_t_' + str(i) + '.fa', 'w') as tf:
                 qf.write('>query\n' + qseq + '\n')
                 tf.write('>target\n' + tseq + '\n')
            process_num += 1

        # lanuch minimap in parallel for this target
        tpath = tmp_path_para + '/' + tname + '_t_{}.fa'
        qpath = tmp_path_para + '/' + tname + '_q_{}.fa'
        apath = tmp_path_para + '/' + tname + '_a_{}.paf'
        process = subprocess.Popen('seq 0 ' + str(process_num - 1) +
                         '|' +  parallel_command + '-I{} -j ' + str(num_thread) + command + tpath + ' ' + 
                         qpath + '''">"''' + apath, shell = True)
        process.wait()
        cs_pattern = re.compile(r'cs:Z:.*?\s')
        for i in range(len(replacement)):
            corrd, qseq = replacement[i]
            _, _, tstart, tend = corrd.get_all() 
            tseq = target_seq_map[tname][tstart:tend]
            res = None
            with open(tmp_path_para + '/' + tname + '_a_' + str(i) + '.paf', 'r') as af:
                res = af.readlines()
            if not res:
                fail_cnt += 1
                print('failed', tname, i, tstart, tend)
                print(qseq)
                print(tseq)
                replacement[i][1] = ''
                continue
            best_alignment = res[0].split()
            qs = int(best_alignment[2])
            qe = int(best_alignment[3])
            ts = int(best_alignment[7])
            te = int(best_alignment[8])
            cs_tag = cs_pattern.findall(res[0])[0].strip()
            qseq = qseq[qs:qe]
            qs = 0
            qe = len(qseq)
            #print('process', i, qseq, cs_tag)
            replace_n = True
            #if replace_n and (qseq.count('N') or qseq.count('n')):
            qseq = parse_cs_tag(qseq, cs_tag, max_length_to_replace)
            qe = len(qseq)
            
            
            while tseq[ts] != qseq[qs] and ts < te and qs < qe:
                ts += 1
                qs += 1
            while tseq[te - 1] != qseq[qe - 1] and te >= ts and qe >= qs:
                te -= 1
                qe -= 1
            
            '''
            while ts < te and qs < qe and tseq[ts] != qseq[qs]:
                ts += 1
                qs += 1
            while te >= ts and qe >= qs and tseq[te - 1] != qseq[qe - 1]:
                te -= 1
                qe -= 1
            '''

            if ts == te or qs == qe:
                replacement[i][1] = ''
                continue
            replacement[i][0].set_all(-1, -1,
                                      tstart + ts, tend - (len(tseq) - te))
            repl_length += tend - (len(tseq) - te) - (tstart + ts)
            replacement[i][1] = qseq[qs:qe]
        shutil.rmtree(tmp_path_para)
                 
    print(repl_length, 'bases has been replaced')
    print('failed', fail_cnt)
    # remove all contents and dirs
    shutil.rmtree(tmp_path)

def out_seq(tname, tseq, fout):
    fout.write('>' + tname + ' len=' + str(len(tseq)) + '\n')
    cnt = 0
    while cnt < len(tseq):
        fout.write(tseq[cnt:cnt + 100] + '\n')
        cnt += 100
    '''    
    for b in tseq:
        fout.write(b)
        cnt += 100
        if cnt == 100:
            fout.write('\n')
    if len(tseq) % 100 != 0:
        fout.write('\n')
    '''
    
def out_new_fasta(replacements, target_seq_map, fasta_out, stat=False):
    length_list = []
    with open(fasta_out, 'w') as fout:
        for tname, tseq in target_seq_map.items():
            if tname in replacements:
                replacement = replacements[tname]
                tseq_new = ''
                pos = 0
                for corrd, seq in replacement:
                    if not seq:
                        continue
                    _, _, tstart, tend = corrd.get_all()
                    tseq_new += tseq[pos:tstart] + seq
                    pos = tend
                tseq_new += tseq[pos:]
                tseq = tseq_new
                if stat:
                    length_list.append((len(tseq), tname))
            out_seq(tname, tseq, fout)
        if stat:
            print('seq stats')
            length_list.sort()
            for slen, sname in length_list:
                print(sname, slen)

def out_new_fasta_debug(replacements, target_seq_map, fasta_out):
    rep_out = fasta_out + '.re.list'
    no_rep_out = fasta_out + '.nore.list'
    rep_length = 0
    norep_length = 0
    with open(rep_out, 'w') as reout, open(no_rep_out, 'w') as noreout:
        for tname, tseq in target_seq_map.items():
            i = 0
            if tname in replacements:
                replacement = replacements[tname]
                pos = 0 # pos on original target
                pos1 = 0 # pos on new target
                for corrd, seq in replacement:
                    if not seq:
                        continue
                    _, _, tstart, tend = corrd.get_all()
                    if pos < tstart:
                        #out_seq(tname + '_' + str(i), tseq[pos:tstart], noreout)
                        length = tstart - pos
                        noreout.write(tname + ':' + str(pos1 + 1) + '-' + str(pos1 + length) + '\n')
                        norep_length += length
                        i += 1
                        pos1 += length
                    #out_seq(tname + '_' + str(i), seq, reout)
                    length = len(seq)
                    reout.write(tname + ':' + str(pos1 + 1) + '-' + str(pos1 + length) + '\n')
                    pos1 += length
                    rep_length += length
                    i += 1
                    pos = tend
                if pos < len(tseq):
                    #out_seq(tname + '_' + str(i), tseq[pos:], noreout)
                    length = len(tseq) - pos
                    noreout.write(tname + ':' + str(pos1 + 1) + '-' + str(pos1 + length) + '\n')
                    norep_length += length
            else:
                #out_seq(tname, tseq, noreout)
                noreout.write(tname + '\n')
                norep_length += len(tseq)
    stat_out = fasta_out + '.stat'
    with open(stat_out, 'w') as sout:
        sout.write('rep_size ' + str(rep_length) + '\n')
        sout.write('unrep_size ' + str(norep_length) + '\n')

def replace_seq(query_fasta_in, target_fasta_in, eval_in, new_target_fasta_out,
                tmp_path, max_length_to_replace, num_thread):
    stat = False
    # create seq mapping to seq name
    query_seq_map = create_seq_map(query_fasta_in)
    target_seq_map = create_seq_map(target_fasta_in, stat)
    
    # debug purpose
    '''
    qseq = query_seq_map['1611']
    tseq = target_seq_map['ctg126']
    print('>1611')
    print(rc(qseq[17617406:17617406 + 132304]))
    print('>ctg126')
    print(tseq[89786:220625])
    exit(0)
    '''
    
    # collect replacements for each target contig
    replacements = parse_eval_file(eval_in, query_seq_map, target_seq_map)
    '''
    for tname, replacement in replacements.items():
        print(tname)
        for corrd, seq in replacement:
            _, _, tstart, tend = corrd.get_all() 
            print(tstart, tend)
            
    exit(1)
    '''

    # rearrage replacements for each target contig
    rearrange_replacements(replacements, target_seq_map, tmp_path,
                           max_length_to_replace, num_thread)

    # output new target fasta file
    out_new_fasta(replacements, target_seq_map, new_target_fasta_out, stat)
    debug = True
    if debug:
        out_new_fasta_debug(replacements, target_seq_map, new_target_fasta_out)
        
    print('successfully finished')


def main():
    # using argparse to parse the command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("query_fasta_in",
                        help="query fasta file", type = str)
    parser.add_argument("target_fasta_in",
                        help="target fasta file", type = str)
    parser.add_argument("eval_in",
                        help="evaluation file", type = str)
    parser.add_argument("new_target_fasta_out",
                        help="new target fasta file", type = str)
    parser.add_argument("--tmp_path", help="new target fasta file",
                        metavar='tmp_path', type = str,
                        default='./tmp/spikeliu/alignment_para')
    parser.add_argument("--max-length-to-replace",
                        help="maximum length to replace(default: 50)",
                        metavar='max_length_to_replace', default=50, type=int)
    parser.add_argument("--num-thread",
                        help="number of threads to use (default: 1)",
                        metavar='num_thread', default=1, type=int)
    args = parser.parse_args()

    replace_seq(args.query_fasta_in, args.target_fasta_in, args.eval_in,
                args.new_target_fasta_out, args.tmp_path,
                args.max_length_to_replace, args.num_thread)


if __name__ == '__main__': 
    main()
