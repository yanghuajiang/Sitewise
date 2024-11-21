import sys, os, subprocess, click, _collections, re, numpy as np, json

def readFasta(fasta, headOnly=False) :
    sequence = _collections.OrderedDict()
    with open(fasta, 'rt') as fin :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                sequence[name] = []
            elif len(line) > 0 and not line.startswith('#') and not headOnly :
                sequence[name].extend(line.strip().split())
    for s in sequence :
        sequence[s] = (''.join(sequence[s])).upper()
    return sequence

def readFastq(fastq) :
    sequence, qual = _collections.OrderedDict(), _collections.OrderedDict()
    with open(fastq, 'rt') as fin :
        line = fin.readline()
        if not line.startswith('@') :
            sequence = readFasta(fastq)
            return sequence, _collections.OrderedDict( [n, re.sub(r'[^!]', 'I', re.sub(r'[^ACGTacgt]', '!', s))] for n, s in sequence.items() )
    with open(fastq, 'rt') as fin :
        for lineId, line in enumerate(fin) :
            if lineId % 4 == 0 :
                name = line[1:].strip().split()[0]
                sequence[name] = []
                qual[name] = []
            elif lineId % 4 == 1 :
                sequence[name].extend(line.strip().split())
            elif lineId % 4 == 3 :
                qual[name].extend(line.strip().split())
    for s in sequence :
        sequence[s] = (''.join(sequence[s])).upper()
        qual[s] = ''.join(qual[s])
    return sequence, qual

complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
def rc(seq, missingValue='N') :
    return ''.join([complement.get(s, missingValue) for s in reversed(seq.upper())])


@click.command()
@click.option('-r', '--ref', help='reference genome in fasta format', required=True)
@click.option('-s', '--site', help='type-specific SNPs in format of <seq_name> <site> <SNP>', required=True)
@click.option('-q', '--qry', help='query MAG in fastq format', required=True)
@click.option('-b', '--bam', help='bam file specifying mapping results', required=True)
@click.option('--minimap2', help='default: /users/softwares/bin/minimap2', default='/users/softwares/bin/minimap2')
@click.option('--samtools', help='default: /users/softwares/bin/samtools', default='/users/softwares/bin/samtools')
def get_site_info(ref, site, qry, bam, minimap2, samtools) :
    sites = {}
    with open(site, 'rt') as fin :
        for line in fin :
            p = line.strip().split()
            if p[0] not in sites :
                sites[p[0]] = []
            if len(p[2]) > 4 or p[2][0] == '.' or p[2][-1] == '.' :
                continue
            sites[p[0]].append([int(p[1]), [p[2][0], p[2][-1]], []])

    qry_seq, _ = readFastq(qry)
    ref_seq, _ = readFastq(ref)


    map_cmd = '{0} -k13 -w5 -c -t1 --frag=yes -A1 -B14 -O24,60 -E2,1 -r100 -g1000 -P -N5000 -f1000,5000 -n2 -m50 -s200 -z200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(
        minimap2, ref, qry
    )

    p = subprocess.Popen(map_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    for line in p.stdout :
        if line.startswith('[') :
            continue
        p = line.strip().split('\t')
        if p[5] not in sites :
            continue
        p[9:11] = int(p[9]), 100. * float(p[9])/float(p[10])
        p[1:4] = [int(p[1]), int(p[2])+1, int(p[3])]
        p[6:9] = [int(p[6]), int(p[7])+1, int(p[8])]

        if p[4] == '+' :
            qi, ri, cigar, d = p[2], p[7], p[-1][5:], 1
        else :
            qi, ri, cigar, d = p[3], p[7], p[-1][5:], -1
        xi = 0
        while xi < len(sites[p[5]]) and sites[p[5]][xi][0] < ri :
            xi += 1
        for s, t in re.findall('(\d+)([MDI])', cigar) :
            s = int(s)
            if t != 'I' :
                rj = ri + s
            if t != 'D' :
                qj = qi + s*d
            while xi < len(sites[p[5]]) and sites[p[5]][xi][0] >= ri and sites[p[5]][xi][0] < rj :
                rd = sites[p[5]][xi][0] - ri
                rx = ri + (rd if t != 'I' else 0) - 1
                qx = qi + (rd*d if t != 'D' else 0) - 1
                rseq = ref_seq[p[5]][rx]
                qseq = qry_seq[p[0]][qx] if d > 0 else rc(qry_seq[p[0]][qx])
                rr = sites[p[5]][xi][2]
                if not (len(rr) and (rr[0][4] >= 10*p[9] or rr[0][5] >= p[10]+0.02)) :
                    sites[p[5]][xi][2].append([p[0], qx+d, rseq, qseq, p[9], p[10]])
                xi += 1
            ri, qi = rj, qj

    base_comp = {}
    if bam :
        qry_sites = {}
        for cont, site_info in sites.items() :
            for site, bb, variants in site_info :
                for var in variants :
                    qry_sites[(var[0], var[1])] = (cont, site)

        with open('{0}.sites'.format(bam), 'wt') as fout :
            for (c, s), (r, s0) in sorted(qry_sites.items()) :
                fout.write('{0}\t{1}\n'.format(c, s))
                fout.write('{0}\t{1}\n'.format(c.split('_', 1)[-1], s))
                fout.write('{0}\t{1}\n'.format(r, s))

        p = subprocess.Popen('{0} mpileup -AB -q 0 -Q 0 --positions {1}.sites {1}'.format(samtools, bam).split(),
                             universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        for line in p.stdout :
            p = line.strip().split('\t')
            bases = re.sub('\^.', '', p[4].upper()).replace('$', '')
            bases = ''.join([b[int(n):] for n, b in re.findall('[+-](\d+)([A-Z]+)', '+0' + bases)])
            bases = dict(zip(*np.unique(list(bases), return_counts=True)))
            base_comp[(p[0], int(p[1]))] = {b:int(bases.get(b, 0)) for b in ('A', 'C', 'G', 'T')}

    print('#Ref_seq\tRef_site\tAnc_base\tAnc_depth\tAlt_seq\tAlt_depth\tDetails')
    for cont, site_info in sites.items() :
        for site, bb, variants in site_info :
            # print(site, bb, variants)
            if variants :
                v = variants[0]
                key = (v[0], v[1])
                if key not in base_comp :
                    key = (v[0].split('_', 1)[-1], v[1])
                if key not in base_comp :
                    key = (cont, v[1])

                bases = base_comp.get(key, {'A':0,'C':0,'G':0,'T':0})
            else :
                bases = {'A':0,'C':0,'G':0,'T':0}
            print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(cont, site, bb[0], bases[bb[0]], bb[1], bases[bb[1]], json.dumps(bases, sort_keys=True)))


if __name__ == '__main__' :
    get_site_info()

