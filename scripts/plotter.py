from __future__ import print_function
import gzip
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import ticker
import seaborn as sns
sns.set_style('white')
sns.set_palette("Set1", 4, 0.75)
colors = sns.color_palette()
from array import array

import collections
import numpy as np
import subprocess as sp
import sys

def fmtr(x, p):
      v, suf = "%.3f" % (x / 1000000.), "M"
      # strip trailing 0's and "."
      while v[-1] in '.0' and '.' in v:
          v = v[:-1]
      return v + suf

def log10(arr):
    arr = np.asarray(arr)
    sign = np.sign(arr)
    arr = np.log10(1 + np.abs(arr))
    return arr * sign

def main(argv):
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument("bigly", nargs="+", help="bigly output(s)")
    p.add_argument("-png", help="Optional PNG output filename")
    a = p.parse_args(argv)


    run(a)

def xopen(f):
    return gzip.open(f) if f.endswith(".gz") else open(f)

header = "chrom pos depth refbase mismatches pairs softstarts softends hardstarts hardends insertstarts insertends deletions splitters splitters1 mean_insert1 mean_insert2 weird discordant discchrom discchromentropy gc65 gc257 duplicity65 duplicity257 spl"
header = header.split()

def run(args):
    gc = 'gc257'

    import matplotlib.gridspec as gridspec

    gs = gridspec.GridSpec(len(args.bigly), 1)
    f = plt.figure(figsize=(9, 5 * len(args.bigly)))
    all_axes = []
    b_axes = []
    l_axes = []

    for bi, bigly in enumerate(args.bigly):
        xs = array('I')
        vals = collections.defaultdict(lambda: array('f'))

        gsi = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[bi], hspace=0.0)
        axes = [plt.Subplot(f, k) for k in gsi]
        for ax in axes:
            f.add_subplot(ax)


        for toks in (l.rstrip("\r\n").split("\t") for l in xopen(bigly)):
            assert len(toks) == len(header), (len(toks), len(header))
            d = dict(zip(header, toks))
            if len(xs) > 1 and int(d['pos']) - 1 != xs[-1]:
                xs.append(xs[-1] + 1)
                xs.append(int(d['pos']) - 1)
                for k in vals:
                    vals[k].extend((0, 0))

            xs.append(int(d['pos']))

            vals['depth'].append(int(d['depth']))
            vals['splitters'].append(int(d['splitters1']))
            vals['softs'].append(int(d['softstarts']) + int(d['softends']))
            vals['discordant'].append(int(d['discordant']))
            if int(d['discchrom']) > 2:
                vals['discchrom'].append(int(d['discchrom'])  * (1 - float(d['discchromentropy'])))
            else:
                vals['discchrom'].append(0)
            vals['inserts1'].append(int(d['mean_insert1']))
            vals['inserts2'].append(int(d['mean_insert2']))
            #vals['duplicity65'].append(float(d['duplicity65']))
            #vals['duplicity257'].append(float(d['duplicity257']))

            vals[gc].append(int(d[gc]))

        for k in vals:
            vals[k] = np.asarray(vals[k])
        xs = np.array(xs)

        #vals['duplicity'] = (vals['duplicity257'] + vals['duplicity65'] ) / 2

        ax = axes[0]
        ax.fill_between(xs, 0, vals['depth'], label='depth', facecolor='gray', lw=0, color='gray', alpha=0.5, zorder=-4)
        ymax=np.percentile(vals['depth'], 99)
        ax.set_ylim(ymax=ymax)
        ax.set_xlim(xs[0], xs[-1])
        all_axes.append(ax)

        #ax.step(xs, vals['duplicity'], lw=1, label='repetitive', color='k')
        ax.step(xs, vals['splitters'], lw=1, label='splitters', color='k')
        #ax.set_ylim(ymax=vals['splitters'].max() + 1)


        s = vals['softs']
        # require at least 10% of reads to be soft-clipped
        s[s < (vals['depth']/10.)] = 0
        xsoft = xs[s > 0]
        ysoft = s[s > 0] #/ vals['depth'][s > 0]
        if len(ysoft) > 0:
            axs = ax.twinx()
            axs.grid('off')
            axs.vlines(xsoft, [0] * len(xsoft), -ysoft, color=colors[2], lw=3, zorder=6, label='soft-clips')
            axs.set_ylim(max(-ymax, 2*-ysoft.max()), 0)
            axs.legend(loc='center right')
        ax.legend()


        ax = axes[1]
        vals['depth'] = vals['depth'].astype(float)

        # plot inserts on a different scale
        axi = ax.twinx()
        axi.step(xs, (vals['inserts1']), label='insert-size-left', ls='none', color=colors[0], alpha=0.65)
        axi.step(xs, (vals['inserts2']), label='insert-size-right', ls='none', color=colors[3], alpha=0.65)
        axi.legend(loc='upper right')
        b_axes.append(axi)

        ax.plot(xs, vals['discordant'], color=colors[1], ls='none', label='discordant')
        ax.step(xs, vals['discchrom'], color=colors[2], ls='none', label='discordant chromosome')
        ax.legend(loc='upper left')
        l_axes.append(ax)

        axes[-1].get_xaxis().set_major_formatter(ticker.FuncFormatter(fmtr))


    ym = max(a.get_ylim()[1] for a in all_axes)
    for a in all_axes:
        a.set_ylim(ymax=ym)
    ym = max(a.get_ylim()[1] for a in b_axes)
    for a in b_axes:
        a.set_ylim(ymax=ym)
    ym = max(a.get_ylim()[1] for a in l_axes)
    for a in l_axes:
        a.set_ylim(ymax=ym)

    plt.subplots_adjust(wspace=0.05, top=0.95, bottom=0.05, right=0.94, left=0.04, hspace=0.11)
    #plt.tight_layout()
    #print(str(args.png))
    if (args.png != None):
        plt.savefig(args.png)
    else:
	#requires X11/GUI
        plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])
