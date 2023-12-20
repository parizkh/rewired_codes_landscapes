import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import logomaker as lm


import gpmap.src.plot.mpl as mplot
import gpmap.src.plot.ds as dplot

from gpmap.src.utils import read_edges
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from gpmap.src.genotypes import get_genotypes_from_region

def rearrange_spines(axes):
    axes.set(aspect='equal', xlabel='', ylabel='')
    axes.spines['left'].set(position=('data', 0), zorder=-1, lw=0.5)
    axes.spines['bottom'].set(position=('data', 0), zorder=-1, lw=0.5)
    # axes.grid(alpha=0.2)
    axes.plot((1), (0), ls="", marker=">", ms=5, color="k",
              transform=axes.get_yaxis_transform(), clip_on=False)
    axes.plot((0), (1), ls="", marker="^", ms=5, color="k",
              transform=axes.get_xaxis_transform(), clip_on=False)
    axes.set(xticks=[-1.0, -0.5, 0.0, 0.5, 1.0],
             yticks=[-1.0, -0.5, 0.0, 0.5, 1.0, 1.5],
             ylim=(None, 1.6), xlim=(None, 1.4))
    
    axes.annotate('Diffusion axis 2', xy=(0.95, 0.07),
                  xycoords=('axes fraction', 'data'), fontsize=8,
                  ha='left', va='center')
    axes.annotate('Diffusion axis 1', xy=(0.075,  0.97),
                  xycoords=('data', 'axes fraction'), fontsize=8,
                  ha='left', va='bottom')
    sns.despine(ax=axes)


if __name__ == '__main__':
    mpl.rcParams['xtick.labelsize'] = 7
    mpl.rcParams['ytick.labelsize'] = 7
    mpl.rcParams['axes.labelsize'] = 8
    mpl.rcParams['axes.titlesize'] = 9

    wt = 'VDGV'
    suffixes = ['FA', 'FG', 'FT', 'FC', 'FV',
                'LA', 'LG', 'LT', 'LC', 'LV',
                'CC', 'CA', 'CG',
                'SA', 'AA', 'AC', 'AG',
                'GA', 'GC', 'GV', 'GL', 'GM', 'GF', 'GT', 'GG', 'GV']
    
    fpath = 'output/protein.nodes.pq'
    print('Loading data from {}'.format(fpath))
    ndf = pd.read_parquet(fpath)
    edf = read_edges('output/protein.edges.npz')

    print('Creating plot for landscape')
    #fig, axes = plt.subplots(1, 1, figsize=(5, 5))
    
    print('\tRendering edges')
    dsg = dplot.plot_edges(ndf, edges_df=edf, x='2', y='1', resolution=2000)
    fig = dplot.dsg_to_fig(dsg)
    axes = fig.axes[0]
    fig.set_size_inches(8, 6)

    cax = inset_axes(axes, width="3%", height="40%",
                         bbox_to_anchor=(-0.05, 0, 1, 0.95),
                         bbox_transform=axes.transAxes, borderpad=0)
    
    print('\tRendering nodes')
    mplot.plot_visualization(axes, ndf, #edges_df=edf,
                             nodes_size=1,
                             nodes_cbar_axes=cax,
                             nodes_cmap_label='Fitness',
                             edges_alpha=0.005,
                             x='2', y='1',
                             nodes_color='function', nodes_cmap='viridis',
                             sort_ascending=True)
    rearrange_spines(axes)
    axes.set(title='')
    
    
    print('Creating boxes around main functional regions')
    boxes = [[(0.5, 1.3), (-0.6, -0.1), 'left'],
             [(-0.4, 0.1), (0.75, 1.48), 'bottom'],
             [(-0.95, -0.05), (-0.8, -0.45), 'bottom']]
    
    for i, (xlims, ylims, pos) in enumerate(boxes):
        mplot.plot_genotypes_box(axes, xlims, ylims, title='',
                                 title_pos=pos, fontsize=8, lw=0.2)
    
    print('Creating logos for selected windows')
    peak_seqs = [get_genotypes_from_region(ndf, min_values={'1': xs[0], '2': ys[0]},
                                           max_values={'1': xs[1], '2': ys[1]})
                 for xs, ys, _ in boxes]
    positions_labels = ['39', '40', '41', '54']
    ndf['region'] = 'None'
    fig2, subplots = mplot.init_fig(3, 1, colsize=2.75, rowsize=1.5)
    for i, (seqs, axes2) in enumerate(zip(peak_seqs, subplots)):
        m = lm.alignment_to_matrix(seqs.values, to_type='probability', pseudocount=0)
        ndf.loc[seqs.values, 'region'] = 'Region {}'.format(i+1)
        m.index = np.arange(m.shape[0])
        logo = lm.Logo(m, ax=axes2, color_scheme='chemistry', vpad=0.05)
        axes2.set(ylabel='Probability', xlabel='Position',
                         xticks=np.arange(m.shape[0]),
                         xticklabels=positions_labels, title='Region {}'.format(i+1))
    fig2.tight_layout()
    print('Rendering and saving logos plot')
    fig2.savefig('plots/protein_landscape.logos.svg', dpi=300)

    print('Rendering and saving visualization plot')
    fig.savefig('plots/protein_landscape.png', dpi=300)

    # Add wt and anchor labels to orient ourselves
    print('\tCreating figure with genotype labels')
    ndf['suffix'] = [x[-2:] for x in ndf.index.values]
    ndf2 = ndf.loc[np.logical_or(np.isin(ndf['suffix'], suffixes), ndf.index.values == wt), :].drop_duplicates(subset='suffix')
    for a, b, label in zip(ndf2['2'], ndf2['1'], ndf2['suffix']):
        axes.text(a, b, label, fontsize=7)
    fig.savefig('plots/protein_landscape.annotated.png', dpi=300)

    print('Writting region classification')
    ndf.to_parquet(fpath)
    
