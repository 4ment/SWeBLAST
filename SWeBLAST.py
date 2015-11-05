import re
import os
import sys
from time import sleep
from argparse import ArgumentParser
import urllib
from subprocess import call
from os.path import expanduser
import ConfigParser

import requests

from SWeBLAST_parser import parse

parser = ArgumentParser(description='Cut sequences in a series of windows and BLAST them using BLAST+.')
parser.add_argument('-i', '--input', metavar='FILE', required=True, help='an input sequence file in FASTA format.')
parser.add_argument('-o', '--output', metavar='FOLDER', required=True, help='an output folder.')
parser.add_argument('-p', '--program', default='blastn',
                    choices=['blastn', 'blastp', 'blastx', 'megablast', 'rpsblast'], help='a blast program.')
parser.add_argument('-d', '--database', default='nr', help='a blast database.')
parser.add_argument('-w', '--window', default='100', type=int, help='window length.')
parser.add_argument('-s', '--step', default='50', type=int, help='step length.')
parser.add_argument('-t', '--threshold', type=float, help='hits with an e-value greater than'
                                                          ' the threshold will be discarded.')
parser.add_argument('-l', '--local', action='store_true', help='use a local database.')

parser.add_argument('-W', '--WWW', action='store_true', help='use old BLAST API instead of BLAST+.')

parser.add_argument('-r', '--parse', action='store_true', help='call SWeBLAST_parser.')
parser.add_argument('-e', '--expect', type=float, help='expect.')
parser.add_argument('-z', '--entrez',  help='ENTREZ query.')

URLBASE = 'http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'

PUT = {'CMD': 'Put',
       'DATABASE': 'nr',
       'PROGRAM': 'blastn'}

GET = {'CMD': 'Get',
       'FORMAT_OBJECT': 'SearchInfo'}

RETRIEVAL = {'CMD': 'Get',
             'FORMAT_TYPE': 'Text'}

PROGRAMS = {'blastp': 'blastp',
            'blastn': 'blastn',
            'blastx': 'blastx',
            'megablast': 'blastn&MEGABLAST=on',
            'rpsblast': 'blastp&SERVICE=rpsblast'}


def get_info(content):
    info = dict()
    for m in re.finditer("<!--.*?QBlastInfoBegin(.*?)QBlastInfoEnd.*?-->", content, re.DOTALL):
        content = m.group(1)
        for mm in re.finditer("\W+(.*?)\W*=\W*(.*)", content):
            info[mm.group(1)] = mm.group(2)
    return info


def lookup_www(query, output_file):
    # build the request
    header = PUT
    header['QUERY'] = urllib.quote_plus(query)

    req = requests.post(URLBASE, params=header)

    if req.status_code == requests.codes.ok:
        info = get_info(req.text)

        # wait for search to complete
        print('Waiting ' + info['RTOE'] + ' seconds RID=' + info['RID'])
        sleep(float(info['RTOE']))

        # poll for results
        while True:
            header = GET
            header['RID'] = info['RID']

            req = requests.post(URLBASE, params=header)

            if req.status_code == requests.codes.ok:
                match = re.search(r'\s+Status=(\w+)', req.text)
                status = match.group(1)

                if status == 'WAITING':
                    match = re.search(r'<p class="WAITING">This page will be automatically updated'
                                      r' in <b>(\d+)</b> seconds</p>', req.text)
                    rtoe = 20
                    if match:
                        rtoe = max(rtoe, int(match.group(1)))

                    print('Waiting ' + str(rtoe) + ' seconds')
                    sleep(rtoe)
                    continue
                elif status == 'READY':
                    if re.search(r'\s+ThereAreHits=yes', req.text) is not None:
                        print('Search complete, retrieving results...')
                        header = RETRIEVAL
                        header['RID'] = info['RID']
                        req = requests.post(URLBASE, params=header)

                        if req.status_code == requests.codes.ok:
                            output_file.write(req.text + '\n')
                        else:
                            sys.stderr.write('Cannot retrieve results.\n' + req.text + '\n')
                    else:
                        print('Search complete, no hits found.')
                elif status == 'FAILED' or status == 'UNKNOWN':
                    sys.stderr.write('Search result ' + status + ' for ' + info['RID'] + '.\n')
                else:
                    sys.stderr.write('Unexpected error during polling\n' + req.text + '.\n')
                break

            # end if is_success
            else:
                sys.stderr.write('Cannot retrieve results during polling\n' + req.text + '\n')

        # end poll while loop
    else:
        sys.stderr.write('Cannot submit query.\nBLAST answer:\n\n' + req.text + '\n')
        if re.search(r'Request-URI\sToo\sLarge', req.text):
            sys.stderr.write('\nSWeBLAST: Use BLAST+ instead.\n\n')


def read_fasta(input_path):
    header = ''
    seq = ''
    sequences = {}
    with open(input_path) as fp:
        for line in fp:
            line = line.rstrip('\n').rstrip('\r')

            if line.startswith('>') and len(seq) > 0:
                sequences[header] = seq
                header = line
                seq = ''
            elif line.startswith('>'):
                header = line
            else:
                seq += line

    sequences[header] = seq
    return sequences


def get_config():
    home = expanduser("~")

    if os.path.exists('.ncbirc'):
        ncbi = '.ncbirc'
    elif os.path.exists('ncbi.ini'):
        ncbi = 'ncbi.ini'
    elif os.path.exists(os.path.join(home, '.ncbirc')):
        ncbi = os.path.join(home, '.ncbirc')
    elif os.path.exists(os.path.join(home, 'ncbi.ini')):
        ncbi = os.path.join(home, 'ncbi.ini')
    else:
        raise NameError('Could not find the NCBI configuration file')

    config = ConfigParser.ConfigParser()
    config.read(ncbi)

    return config


def main():
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    output_path = os.path.join(args.output, 'blast.txt')
    output_file = open(output_path, 'w')

    sequence_path = os.path.join(args.output, 'seq.fa')
    sequence_file = open(sequence_path, 'w')

    cmd = ''

    if args.WWW:
        if args.entrez is not None:
            PUT['ENTREZ_QUERY'] = args.entrez
        if args.database is not None:
            PUT['DATABASE'] = args.database
        PUT['PROGRAM'] = PROGRAMS[args.program]
    else:
        print cmd
        cmd = args.program + ' -db ' + args.database + ' -out ' + output_path + ' -query ' + sequence_path
        if args.local:
            config = get_config()
            print('Using local BLAST config file: ' +os.path.join(config.get('BLAST', 'BLASTDB'), args.database))
        else:
            cmd += ' -remote'

    window = args.window
    step = args.step

    sequences = read_fasta(args.input)

    for name in sequences.keys():
        if args.step == 0 and window == 0:
            query = name + '\n' + sequences[name]
            lookup_www(query, output_file)
        else:
            i = 0
            while i+window < len(sequences[name]):
                begin = i + 1
                end = i + window
                sequence = sequences[name][i:i+window]
                query = name + '|' + str(begin) + '-' + str(end) + '\n' + sequence
                sequence_file.write(query + '\n')

                if args.WWW:
                    print(query)
                    lookup_www(query, output_file)

                i += step

    sequence_file.close()
    output_file.close()

    if not args.WWW:
        print('CMD: ' + cmd)
        return_code = call(cmd, shell=True)
        if return_code != 0:
            print('The following command failed: ' + cmd)

    if args.parse:
        parse(output_path, os.path.join(args.output, 'blast_STAT'), args.threshold)

main()
