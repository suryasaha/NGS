#!/usr/bin/python

#Modified from Sylvan's libngs script
#https://github.com/sylvainforet/libngs/blob/master/src/scripts/assembly_stats.py

import cStringIO

statsNames = ['name', 'Contigs', 'totLength', 'Longest', 'Shortest', 'Ns',
              'N75', 'N50', 'N25',
              'len75', 'len50', 'len25',
              'count400', 'count1k', 'count5k','count10k','count20k','count40k','count60k','count70k','count80k','count90k',
              'sum1k', 'sum5k', 'sum10k', 'sum20k', 'sum40k', 'sum60k', 'sum70k', 'sum80k', 'sum90k']

def loadAssembly(path):
    handle = open(path)
    seqs = {}
    name = None
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':
            name = line
            seqs[name] = cStringIO.StringIO()
        elif name:
            seqs[name].write(line)
    handle.close()
    for name in seqs:
        tmp = seqs[name].getvalue()
        seqs[name].close()
        seqs[name] = tmp
    sizes = [len(x) for x in seqs.values()]
    sizes.sort()
    return sizes

def getNcounts(path):
    handle = open(path)
    seqs = {}
    name = None
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':
            name = line
            seqs[name] = cStringIO.StringIO()
        elif name:
            seqs[name].write(line)
    handle.close()
    for name in seqs:
        tmp = seqs[name].getvalue()
        seqs[name].close()
        seqs[name] = tmp
    ncounts = [x.lower().count('n') for x in seqs.values()]
    return ncounts

def countLargerThan(sizes, minimum):
    largerSeqs = [x for x in sizes if x >= minimum]
    return len(largerSeqs)

def sumLargerThan(sizes, minimum):
    largerSeqs = [x for x in sizes if x >= minimum]
    return sum(largerSeqs)

def printHeader(handle):
    header = '\t'.join(statsNames)
    handle.write(header + '\n')

def getN25_50_75(sizes):
    total = sum(sizes)
    tmp = 0
    n25 = 0
    n50 = 0
    n75 = 0
    for i in sizes:
        tmp += i
        if tmp >= 0.25 * total and not n75:
            n75 = i
        if tmp >= 0.50 * total and not n50:
            n50 = i
        if tmp >= 0.75 * total and not n25:
            n25 = i
    return n25, n50, n75

def printAssemblyStats(ncounts, sizes, name, handle):
    nContigs = len(sizes)
    stats = {}
    stats['Ns' ] = sum(ncounts)
    stats['name' ] = name
    stats['len75' ] = sizes[nContigs / 4]
    stats['len50' ] = sizes[nContigs / 2]
    stats['len25' ] = sizes[(3 * nContigs) / 4]
    a, b, c = getN25_50_75(sizes)
    stats['N25' ] = a
    stats['N50' ] = b
    stats['N75' ] = c
    stats['Contigs' ] = nContigs
    stats['Longest' ] = max(sizes)
    stats['Shortest' ] = min(sizes)
    stats['count400' ] = countLargerThan(sizes, 400)
    stats['count1k'] = countLargerThan(sizes, 1000)
    stats['count5k'] = countLargerThan(sizes, 5000)
    stats['count10k'] = countLargerThan(sizes, 10000)
    stats['count20k'] = countLargerThan(sizes, 20000)
    stats['count40k'] = countLargerThan(sizes, 40000)
    stats['count60k'] = countLargerThan(sizes, 60000)
    stats['count70k'] = countLargerThan(sizes, 70000)
    stats['count80k'] = countLargerThan(sizes, 80000)
    stats['count90k'] = countLargerThan(sizes, 90000)
    stats['totLength' ] = sumLargerThan(sizes, 0)
    stats['sum1k' ] = sumLargerThan(sizes, 1000)
    stats['sum5k' ] = sumLargerThan(sizes, 5000)
    stats['sum10k'] = sumLargerThan(sizes, 10000)
    stats['sum20k'] = sumLargerThan(sizes, 20000)
    stats['sum40k'] = sumLargerThan(sizes, 40000)
    stats['sum60k'] = sumLargerThan(sizes, 60000)
    stats['sum70k'] = sumLargerThan(sizes, 70000)
    stats['sum80k'] = sumLargerThan(sizes, 80000)
    stats['sum90k'] = sumLargerThan(sizes, 90000)
    handle.write(name + '\t')
    for i in statsNames[1:-1]:
        handle.write('%d\t' % stats[i])
    handle.write('%d\n' % stats[statsNames[-1]])

def main():
    import sys
    import os.path
    if len(sys.argv) < 2:
        print 'Usage: %s FILE' % os.path.basename(sys.argv[0])
        sys.exit(1)
    printHeader(sys.stdout)
    for i in sys.argv[1:]:
        sizes = loadAssembly(i)
        ncounts = getNcounts(i)
        printAssemblyStats(ncounts, sizes, i, sys.stdout)

if __name__ == '__main__':
    main()

###
