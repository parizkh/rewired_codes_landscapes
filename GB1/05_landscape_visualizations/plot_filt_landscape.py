import sys
import numpy as np
import pandas as pd
import holoviews as hv
import seaborn as sns
import matplotlib.pyplot as plt

import gpmap.src.plot as plot

from gpmap.src.utils import read_dataframe, read_edges


def make_plot(ndf, edf, pathdf, x, y, z, code, suffixes, size=5, padding=0.1):
    xlim = ndf[x].min(), ndf[x].max()
    ylim = ndf[y].min(), ndf[y].max()
    dx, dy = xlim[1]- xlim[0], ylim[1] - ylim[0]
    
    # Plot edges
    print('\tCreating edges')
    dsg = plot.plot_edges_datashader(ndf, edges_df=edf, x=x, y=y,
                                         resolution=1500)
    dsg.opts(aspect='square', padding=0.1)
    fig = hv.render(dsg)
    fig.set_size_inches(size, size)
    axes = fig.axes[0]
    
    # Plot nodes
    print('\tnodes')
    plot.plot_nodes(axes, ndf, size=2, vmin=vmin, vmax=vmax, cmap=cmap, color='function', ascending=True,
                    x=x, y=y, sort_by=z, cbar=False, clabel='log(relative binding)')

    # Fine tuning
    ylim = (ylim[0] - padding * dy, ylim[1] + padding * dy)
    xlim = (xlim[0] - padding * dx, xlim[1] + padding * dx)
    axes.set(xlim=xlim, ylim=ylim)

    # Saving figure
    print('\tRendering unannotated figure')
    fig.savefig('plots/{}_{}.{}.filtered.png'.format(code, x, y), dpi=300)
    
    if code == 'Standard':
        print('\tCreate figure with path')
        plot_path(axes, pathdf, x, y, vmin, vmax, cmap)
        fig.savefig('plots/{}_{}.{}.path.png'.format(code, x, y), dpi=300)
    
    # Add wt and anchor labels to orient ourselves
    ndf['suffix'] = [x[-2:] for x in ndf['protein']]
    ndf2 = ndf.loc[np.isin(ndf['suffix'], suffixes), :].drop_duplicates(subset='suffix')
    for a, b, label in zip(ndf2[x], ndf2[y], ndf2['suffix']):
        axes.text(a, b, '--' + label, fontsize=7)
        
    ndf2 = ndf.loc[ndf['protein'] == wt, :].drop_duplicates(subset='suffix')
    plot.plot_nodes(axes, ndf2, x=x, y=y, size=20, 
                    vmin=vmin, vmax=vmax, cmap=cmap, color='function', cbar=False)
    for a, b, label in zip(ndf2[x], ndf2[y], ndf2['protein']):
        axes.text(a, b, label, fontsize=8)
        
    # Saving figure
    print('\tRender figure with genotype annotations')
    fig.savefig('plots/{}_{}.{}.annotated.png'.format(code, x, y), dpi=300)


def read_landscape(code, seqs):
    fpath = 'output/{}.nodes.pq'.format(code)
    df = read_dataframe(fpath)
    wtdf = df.loc[df['protein'] == wt, :]
    pathdf = df.loc[seqs, :]
    
    fpath = 'output/{}.filtered.nodes.pq'.format(code)
    ndf = read_dataframe(fpath)
    
    ndf = pd.concat([ndf, wtdf])
    ndf.to_parquet(fpath)
    
    fpath = 'output/{}.filtered.edges.npz'.format(code)
    edf = read_edges(fpath)
    return(ndf, edf, pathdf)


def get_path_seqs():
    background = 'GUCGAU'
    seqs = ['GGUCUU', 'GGUGUU', 'GGUGCU', 'GCUGCU',
            'UCUGCU', 'UGUGCU', 'UGUGGU', 'UUUGGU',
            'CUUGGU', 'CUUGCU', 'CUUACU']
    l = len(seqs)
    genotypes = [background + seq for seq in seqs]
    genotypes = [x.replace('U', 'T') for x in genotypes]
    return(genotypes)


def plot_path(axes, ndf, x, y, vmin, vmax, cmap):
    l = ndf.shape[0]
    edf = pd.DataFrame({'i': np.arange(l-1), 'j': np.arange(1, l)})
    plot.plot_nodes(axes, ndf, x, y, size=40, cmap=cmap, cbar=False,
                    zorder=4, lw=1.5, vmax=vmax, vmin=vmin)
    plot.plot_edges(axes, ndf, edf, x, y, width=1.5, alpha=1, zorder=3,
                    color='black')


if __name__ == '__main__':
    vmin, vmax = (-8.078194290399551, 2.297213241457939)
    cmap = 'viridis'
    wt = 'VDGV'
    suffixes = ['FA', 'FG', 'FT', 'FC', 'FV',
                'LA', 'LG', 'LT', 'LC', 'LV',
                'CC', 'CA', 'CG',
                'SA', 'AA', 'AC', 'AG',
                'GA', 'GC', 'GV', 'GL', 'GM', 'GF', 'GT', 'GG', 'GV']
    seqs = get_path_seqs()

    code, x, y, z = 'Standard', '2', '1', '3'
    #code = sys.argv[1]
    
    print('Reading filtered landscape under code {}'.format(code))
    ndf, edf, pathdf = read_landscape(code, seqs)
    make_plot(ndf, edf, pathdf, x, y, z, code, suffixes, size=5, padding=0.1)


