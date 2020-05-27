import pandas as pd
from optparse import OptionParser
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import subprocess

def getOptions():
    parser = OptionParser()

    parser.add_option("--f", dest = "infile",
        help = "File of distances from closest CAGE peak", metavar = "FILE", type = str)
    parser.add_option("--xmax", dest = "xmax", type = int,
                      help = "Max x value for plots.", default = 5000)
    parser.add_option("--ymax", dest = "ymax", type = int,
                      help = "Max y value for plots.", default = 1000000)
    parser.add_option("--o", dest = "outprefix", help = "Prefix for outfile",
        metavar = "FILE", type = str)

    (options, args) = parser.parse_args()
    return options

def plot_histogram(data, xvar, label, xmax, ymax, fname):
    x = pd.Series(data[xvar], name=label)
    plt.xlim(-1*xmax,xmax)
    plt.ylim(0,ymax)
    ax = sns.distplot(x, kde = False, color = "dodgerblue",
                  bins = np.arange(-1*xmax -1, xmax + 1, 5))

    med = round(np.median(x), 1)

    style = dict(size=12, color='black')
    plt.axvline(med, linestyle = '--', color = 'lightgrey')

    ax.text(med + 50, ymax*7/8, "Median: " + str(med) + " bp", **style)

    plt.tight_layout()
    plt.savefig(fname, dpi = 600, bbox_inches='tight')
    plt.close()


def main():
    options = getOptions()

    data = pd.read_csv(options.infile, sep = '\t', header = None)
    data = data.rename(columns={ data.columns[-1]: "dist" })
    print(data)
    
    fname = options.outprefix + "_closest_CAGE.png"
    label = "Distance from read start to closest CAGE peak (bp)"
    plot_histogram(data, "dist", label, options.xmax, options.ymax, fname)
   
if __name__ == '__main__':
    main() 
