#!/usr/bin/python
r""" 
    Plot the Orion proto transition matrix for a certain level
"""

import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from operator import itemgetter


def get_args():
    r""" Parse arguments for the program    
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("transfile", help="File containing transition probabilities for a given transition level.")
    parser.add_argument('-ns', '--noshow', action='store_true', default = False, help = 'Do not display the figure.')
    parser.add_argument("-t","--title", type=str, default=None, help="Chart title.")
    parser.add_argument("-o","--ofile", type=str, default=None, help="Output file for the chart (must have an accepted pyplot output format extension).")
    parser.add_argument("-m","--minp", type=float, default=0.01, help="Minimum transition probability to display text for.")
    
    return parser.parse_args()


if __name__ == '__main__':

    args = get_args()  # command line arguments
    
    data = []
    edata = []
    
    # load data file
    with open(args.transfile,"r") as fh:
        for line in fh:
            p = line.split()
            assert len(p) == 3
            if p[0] == 'S':
                data.append((0, int(p[1])+1, float(p[2])))
            elif p[1] == 'E':
                edata.append((int(p[0])+1, 'E', float(p[2])))
            else:
                data.append((int(p[0])+1, int(p[1])+1, float(p[2])))
    
    # find number of rows/cols
    nrows = max(data, key=itemgetter(0))[0] + 2
    ncols = max(data, key=itemgetter(1))[1] + 2
    nrows = ncols = max(nrows, ncols)
    
    # transfer data to matrix
    mat = np.zeros((nrows, ncols), dtype=np.float)
    for f in data:
            mat[f[0], f[1]] = f[2]
    for f in edata:
            mat[f[0], ncols-1] = f[2]
    
    # create figure
    fig, ax = plt.subplots()
    cm = plt.pcolor(mat, cmap=plt.cm.Blues, edgecolors = 'None', )
    labels = ['S']
    nrows = mat.shape[0]
    labels.extend(range(nrows-2))
    labels.append('E')
    
    plt.xticks(np.arange(0.5, nrows+0.6), labels)
    plt.yticks(np.arange(0.5, ncols+0.6), labels)
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.gca().set_xlim((0, nrows))
    plt.gca().set_ylim((0, ncols))
    
    # add text to cells
    for i in xrange(nrows):
        for j in xrange(ncols):
            if mat[i][j] >= args.minp:
                cl = '0.2' if mat[i][j] <= 0.4 else '0.8'
                ax.text(j+0.5, i+0.5, "%s" % ("%.2f" % mat[i][j])[1:] if mat[i][j] < 1 else "%.1f" % mat[i][j], 
                        va='center', ha='center', fontsize=8, color=cl )    
    
    fig.set_size_inches(7.3, 6)
    if args.title:
        plt.title(args.title, fontsize=24)
    plt.tight_layout()
    if args.ofile:
        plt.savefig(args.ofile, dpi=300)
        print "Saved figure to '%s'." % args.ofile
    if not args.noshow:
        plt.show()
