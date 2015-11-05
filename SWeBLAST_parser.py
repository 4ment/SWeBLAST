import re
import os
from argparse import ArgumentParser
from math import log, sqrt, isnan
import copy


def parse(input_path, output_path, threshold):

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    summary_file = open(os.path.join(output_path, 'summary.csv'), 'w')

    re_program = re.compile(r'^(t?blast[pnx])\s.+', re.IGNORECASE)
    re_query = re.compile(r'^Query=\s+(\S.+)')
    re_hit = re.compile('^>(.+)')
    re_score = re.compile(r'^\s*Score\s+=.+Expect\s*=\s*(\d.+)')
    re_frame = re.compile(r' Frame = (.+)')

    nb = 0
    e = ''
    previous = ''
    name = ''
    ranking = 0
    program = ''

    hits = {}

    with open(input_path) as fp:
        for line in fp:
            line = line.rstrip('\n').rstrip('\r')

            # Get program: blastn, blastp, blastx
            match_hit = re_program.match(line)
            if match_hit:
                program = match_hit.group(1)

            # Get query name
            match_hit = re_query.match(line)
            if match_hit:
                seq_query = match_hit.group(1)
                nb += 1
                ranking = 0
                previous = ""
                for l in fp:
                    if re.match('^Length', l):
                        break
                    seq_query += ' ' + l

                seq_query = seq_query.replace(',', '')
                seq_query = seq_query.rstrip('\n').rstrip('\r')
                summary_file.write(seq_query)
                continue

            # HIT
            match_hit = re_hit.match(line)
            match_frame = re_frame.match(line)

            if re_hit.match(line):
                name = match_hit.group(1)
                for l in fp:
                    if re.match('^Length', l):
                        break
                    name += ' ' + l
                    name = name.rstrip('\n').rstrip('\r')
                name = name.replace(',', '')

                for l in fp:
                    match_hit = re_score.match(l)
                    if match_hit:
                        e = match_hit.group(1)
                        e = e.rstrip('\n').rstrip('\r')
                        break
                summary_file.write(',' + name + '(' + e + ')')

                if previous == '' or previous != e:
                    ranking += 1  # In case of a tie

                if threshold is None or (threshold is not None and e < threshold):
                    if name in hits:
                        if nb >= len(hits[name]['EVALUE']):
                            hits[name]['EVALUE'].extend([-1 for x in xrange(len(hits[name]['EVALUE']), nb+1)])
                            hits[name]['RANK'].extend([-1 for x in xrange(len(hits[name]['RANK']), nb+1)])
                        hits[name]['EVALUE'][nb] = float(e)
                        hits[name]['RANK'][nb] = int(ranking)
                    else:
                        es = rankings = [-1 for x in range(nb+1)]
                        es[nb] = float(e)
                        rankings[nb] = int(ranking)
                        hits[name] = {'EVALUE': es, 'RANK': rankings}

                previous = e

            # Used when no hits found
            elif re.match(r'<b>No significant similarity found.</b>', line):
                summary_file.write(',No hits found')

            # For BLASTX
            elif match_frame:
                # if previously added then it is above the threshold
                if name in hits:
                    hits[name]['FRAME'][nb] = match_frame.group(1)

            # End of one blast
            elif re.match('^  Database', line) or re.match('^\s*Effective search space used', line):
                summary_file.write('\n')

    summary_file.close()

    order_log_e_values(hits)
    order_ranks(hits)

    print_e_values(hits, output_path)
    print_ranks(hits, output_path)
    print_stats(hits, output_path)
    if re.match('blastz', program, re.IGNORECASE):
        print_frames(hits, output_path)


def order_log_e_values(hits):
    for f in hits.keys():
        e_values = copy.copy(hits[f]['EVALUE'])
        while -1 in e_values:
            e_values.remove(-1)
        for i in range(len(e_values)):
            if e_values[i] == 0:
                e_values[i] = 9999
            else:
                e_values[i] = abs((log(e_values[i]) / log(10)))

        hits[f]['EVALUEORDERED'] = sorted(e_values)


def order_ranks(hits):
    for f in hits.keys():
        rank = copy.copy(hits[f]['RANK'])
        while -1 in rank:
            rank.remove(-1)

        hits[f]['RANKORDERED'] = sorted(rank)


def print_stats(hits, output):

    stat_file = open(os.path.join(output, 'stat.csv'), 'w')

    stat_file.write('Sequence,nb,min e-value,max e-value,range e-value,median e-value,mean e-value,std e-value,'
                    'min rank,max rank,range rank,median rank,mean rank,std rank\n')

    for f in hits.keys():
        # e-value
        e_values = hits[f]['EVALUEORDERED']

        e_range = abs(e_values[0] - e_values[len(e_values) - 1])
        e_mean = mean(e_values)
        e_std = std(e_values, e_mean)

        # Rank
        rank = hits[f]['RANKORDERED']

        rank_range = abs(rank[len(rank) - 1] - rank[0])
        rank_mean = mean(rank)
        rank_std = std(rank, rank_mean)

        line = f + ',' + str(len(e_values)) + ',' + str(e_values[len(e_values) - 1]) + ',' + str(e_values[0]) + ','
        line += str(e_range) + ',' + str(median(e_values)) + ',' + str(e_mean) + ',' + str(e_std) + ','
        line += str(rank[0]) + ',' + str(rank[len(rank) - 1]) + ',' + str(rank_range) + ',' + str(median(rank)) + ','
        line += str(rank_mean) + ',' + str(rank_std) + ',\n'
        stat_file.write(line)

    stat_file.close()


def print_e_values(hits, output):
    e_value_file = open(os.path.join(output, 'evalue.csv'), 'w')
    for f in hits.keys():
        e_value_file.write(f + ',')
        for e in hits[f]['EVALUE']:
            if e != - 1:
                if e == 0:
                    e_value_file.write('99999')
                else:
                    e_value_file.write(str(abs(log(e) / log(10))))
            e_value_file.write(',')
        e_value_file.write('\n')

    e_value_file.close()


def print_ranks(hits, output):
    rank_file = open(os.path.join(output, 'rank.csv'), 'w')
    for f in hits.keys():
        rank_file.write(f + ',')
        for r in hits[f]['RANK']:
            if r != -1:
                rank_file.write(str(r))
            rank_file.write(',')
        rank_file.write('\n')
    rank_file.close()


def print_frames(hits, output):
    frame_file = open(os.path.join(output, 'frame.csv'), 'w')
    for f in hits.keys():
        frame_file.write(f + ',')
        for r in hits[f]['FRAME']:
            if r != -1:
                frame_file.write(r)
            frame_file.write(',')
        frame_file.write('\n')
    frame_file.close()


def median(x):
    return x[int(len(x)/2)]


def mean(x):
    return float(sum(x))/len(x) if len(x) > 0 else float('nan')


def std(x, avg):
    var = variance(x, avg)
    return sqrt(var) if not isnan(var) else var


def variance(x, avg):
    length = len(x)
    var = float('nan')
    if length == 1:
        var = 0
    elif length > 1:
        ssq = map(lambda y: (y - avg)**2, x)
        var = sum(ssq) / (length - 1)
    return var


if __name__ == '__main__':
    parser = ArgumentParser(description='Parse a blast file generated from SWeBLAST.pl and create summary files.')
    parser.add_argument('-i', '--input', metavar='FILE', required=True, help='an input file.')
    parser.add_argument('-o', '--output', metavar='FOLDER', required=True, help='an output folder.')
    parser.add_argument('-t', '--threshold', type=float, help='hits with an e-value greater than '
                                                              'the threshold will be discarded.')
    args = parser.parse_args()
    parse(args.input, args.output, args.threshold)
