#!/share/app/python-3.4.3/bin/python3
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "spikeliu"

import argparse
import re

def create_dummy_files(quast_stdout_file,
                       dummy_blat_file,
                       dummy_scaff_file):
    with open(quast_stdout_file, 'r') as quast_in_f, \
         open(dummy_blat_file, 'w') as blat_out_f, \
         open(dummy_scaff_file, 'w') as scaff_out_f:
         
         contig_cnt = scaff_cnt = 0
         prev_strand = ''
         strand_error = error_cnt = 0
         strand_error_v1 = 0
         forward = {}
         backward = {}
         cnt = 0
         inversion_len = 0
         for line in quast_in_f:
             line = line.strip()
             if line.startswith('CONTIG: '):
                 scaff_cnt += 1
                 #scaff_name = '>scaffold' + str(scaff_cnt)
                 match = re.match(r'^CONTIG: (.+?) \(([0-9]+)bp\)$', line)
                 scaff_name = '>' + match.group(1) + ' ' + match.group(2)
                 scaff_out_f.write(scaff_name + '\n')
                 if cnt == 0:
                     continue
                 # print ('nimei ' + str(error_cnt))
                 if error_cnt % 2 == 0:
                     strand_error += error_cnt / 2
                 else:
                     strand_error += error_cnt / 2 + 1
                 error_cnt = 0
                 prev_strand = ''
                 # suppose forward is always more than backward
                 if len(forward) < len(backward):
                     backward, forward = forward, backward
                 if len(forward) != cnt:
                     for key, value in backward.items():
                         #print(key, value)
                         inversion_len += value
                     strand_error_v1 += len(backward)
                 forward = {}
                 backward = {}
                 cnt = 0
             elif line.startswith('Real Alignment') or \
                  line.startswith('One align') or \
                  line.startswith('Alignment'):
                 contig_cnt += 1
                 cnt += 1
                 qname = str(contig_cnt)
                 content = re.split(':|\|', line)
                 # start and end pos on ref
                 tstart, tend = (int(x) for x in content[1].strip().split())
                 tstart -= 1
                 qstart, qend = (int(x) for x in content[2].strip().split())
                 # reverse complement contig
                 strand = '+'
                 if qstart > qend:
                     strand = '-'
                     qend, qstart = qstart, qend
                 offset = qstart - 1
                 qstart = 0
                 qend -= offset
                 match, qlen = (int(x) for x in content[3].strip().split())

                 # record inversion occurence
                 if strand == '+':
                     forward[qname] = qlen 
                 else:
                     backward[qname] = qlen

                 tname = content[5].strip().split()[0]
                 # create content for blat output
                 content = [str(qlen), 'NULL', 'NULL', 'NULL', 'NULL', '0',
                            'NULL', '0', '+', qname, str(qlen), str(qstart),
                            str(qend), tname, 'NULL', str(tstart), str(tend),
                            'NULL', 'NULL', 'NULL']
                 blat_out_f.write('\t'.join(content) + '\n')
                 # create content for scaff output
                 content = [qname, str(offset), strand, str(qlen)]
                 scaff_out_f.write(' '.join(content) + '\n')

                 # check if strand is wrong
                 if prev_strand:
                     if prev_strand != strand:
                         error_cnt += 1
                 # record previous strand
                 prev_strand = strand
         # deal with last scaff
         if cnt > 0:
             if error_cnt % 2 == 0:
                 strand_error += error_cnt / 2
             else:
                 strand_error += error_cnt / 2 + 1

             error_cnt = 0
             prev_strand = ''
             # suppose forward is always more than backward
             if len(forward) < len(backward):
                 backward, forward = forward, backward
             if len(forward) != cnt:
                 for key, value in backward.items():
                     #print(key, value)
                     inversion_len += value
                 strand_error_v1 += len(backward)
             forward = {}
             backward = {}
             cnt = 0
         #print('strand error: ', strand_error)
         #print('strand error v1: ', strand_error_v1)
         #print('inversion length: ', inversion_len)

def parse_contig_distance_on_ref_in_quast(quast_file, out_file):
    lines = []
    with open(quast_file, 'r') as q_f:
        for line in q_f:
            lines.append(line.strip())
    
    with open(out_file, 'w') as o_f:
        for i in range(len(lines)):
            line = lines[i]
            if line.startswith('Incorrectly estimated size'):
                 for j in reversed(range(i)):
                     if lines[j].startswith('Real Alignment') or \
                         lines[j].startswith('One align') or \
                         lines[j].startswith('Alignment'):
                         break
                 content = re.split(':|\|', lines[j])
                 # start and end pos on ref
                 tstart_lhs, tend_lhs = (int(x) for x in content[1].strip().split())
                 qstart_lhs, qend_lhs = (int(x) for x in content[2].strip().split())
                 for j in range(i + 1, len(lines)):
                     if lines[j].startswith('Real Alignment') or \
                         lines[j].startswith('One align') or \
                         lines[j].startswith('Alignment'):
                         break
                 content = re.split(':|\|', lines[j])
                 # start and end pos on ref
                 tstart_rhs, tend_rhs = (int(x) for x in content[1].strip().split())
                 qstart_rhs, qend_rhs = (int(x) for x in content[2].strip().split())
                 if qstart_lhs < qend_lhs and qstart_rhs < qend_rhs:
                     o_f.write(str(tstart_rhs - tend_lhs) + '\n')
                 elif qstart_lhs > qend_lhs and qstart_rhs > qend_rhs:
                     o_f.write(str(tstart_lhs - tend_rhs) + '\n')
                 else:
                     print('that should not happen.')
                     print(lines[i - 1])
                     print(lines[i])
                     print(lines[i + 1])

def parse_contig_in_quast(quast_file, out_file, min_length):
    lines = []
    with open(quast_file, 'r') as q_f:
        for line in q_f:
            lines.append(line.strip())
    with open(out_file, 'w') as o_f:
        for i in range(len(lines)):
            line = lines[i]
            if line.startswith('Real Alignment') or \
                line.startswith('One align') or \
                line.startswith('Alignment'):
                content = re.split(':|\|', line)
                # start and end pos on ref
                tstart, tend = (int(x) for x in content[1].strip().split())
                tstart -= 1
                if tend - tstart < min_length:
                    continue
                tname = content[5].strip().split()[0]
                o_f.write(str(tname) + '\t' + str(tstart) + '\t' + str(tend) + '\n')
                     
def main():
    # using argparse to parse the command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("quast_file",
                        help="quast input, usuall is"
                        "<quast_out_root>/contigs_reports/"
                        "contig_reports_xxx.stdout", type = str)
    parser.add_argument("scaff_file", help="scaffold output", type = str)
    parser.add_argument("blat_file", help="blat output", type = str)
    args = parser.parse_args()
    create_dummy_files(args.quast_file,
                       args.blat_file,
                       args.scaff_file)

if __name__ == '__main__': 
    main()
