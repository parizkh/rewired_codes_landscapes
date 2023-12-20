import sys
import numpy as np
import pandas as pd
import holoviews as hv
import matplotlib.cm as cm

import gpmap.src.plot.mpl as mplot
import gpmap.src.plot.ds as dplot

from gpmap.src.utils import read_dataframe, read_edges

import matplotlib.colors as mc
import colorsys

def adjust_lightness(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def make_plot(ndf, edf, x, y, code, suffixes, size=5, padding=0.1, pathdf=None):
    showlegend = code == 'Standard' and x == '2'
    cmap = cm.get_cmap('binary')
    palette = {'None': '#7E7E7E',     # cmap(0.2),
               'Region 1': '#648FE4', # 'darkred',
               'Region 2': '#F2D022', #'orange',
               'Region 3': '#5B44AF', #'purple'
               }
    
    xlim = ndf[x].min(), ndf[x].max()
    ylim = ndf[y].min(), ndf[y].max()
    dx, dy = xlim[1]- xlim[0], ylim[1] - ylim[0]
    
    # Plot edges
    print('\tCreating edges')
    dsg = dplot.plot_edges(ndf, edges_df=edf, x=x, y=y, resolution=800)
    dsg.opts(aspect='square', padding=0.1)
    fig = hv.render(dsg)
    fig.set_size_inches(size, size)
    axes = fig.axes[0]
    
    # Plot nodes
    print('\tCreating nodes')
    ndf['color'] = [adjust_lightness(palette[x], amount=np.random.normal(1, 0.05)) for x in ndf['region']]

    ndf_aux = ndf.sort_values('function', ascending=True)
    axes.scatter(ndf_aux[x], ndf_aux[y], s=2, c=ndf_aux['color'])
    #mplot.plot_nodes(axes, ndf, x=x, y=y, size=2, 
    #                 color='color', palette=palette, legend=False,
    #                 sort_by='function', sort_ascending=True)
    if showlegend:
        colors_dict = {k: v for k, v in palette.items() if k != 'None'}
        print(colors_dict)
        mplot.create_patches_legend(axes, colors_dict, loc=4, fontsize=8)

    # Fine tuning
    ylim = (ylim[0] - padding * dy, ylim[1] + padding * dy)
    xlim = (xlim[0] - padding * dx, xlim[1] + padding * dx)
    axes.set(xlim=xlim, ylim=ylim,
             xlabel='Diffusion axis {}'.format(x),
             ylabel='Diffusion axis {}'.format(y))

    # Saving figure
    print('\tRendering unannotated figure')
    fig.savefig('plots/{}_{}.{}.filtered.png'.format(code, x, y), dpi=300)
    
    # Creating additional figure with specified path
    if pathdf is not None:
        print('\tCreate figure with path')
        plot_path(axes, pathdf, x, y, vmin, vmax, cmap)
        fig.savefig('plots/{}_{}.{}.path.png'.format(code, x, y), dpi=300)
    
    # Add wt and anchor labels to orient ourselves
    print('\tCreating figure with genotype labels')
    ndf['suffix'] = [x[-2:] for x in ndf['protein']]
    ndf2 = ndf.loc[np.isin(ndf['suffix'], suffixes), :].drop_duplicates(subset='suffix')
    for a, b, label in zip(ndf2[x], ndf2[y], ndf2['suffix']):
        axes.text(a, b, '--' + label, fontsize=7)
        
    ndf2 = ndf.loc[ndf['protein'] == wt, :].drop_duplicates(subset='suffix')
    mplot.plot_nodes(axes, ndf2, x=x, y=y, size=20, 
                     vmin=vmin, vmax=vmax, cmap=cmap, color='function', cbar=False)
    for a, b, label in zip(ndf2[x], ndf2[y], ndf2['protein']):
        axes.text(a, b, label, fontsize=8)
        
    # Saving figure
    print('\tRender figure with genotype annotations')
    fig.savefig('plots/{}_{}.{}.annotated.png'.format(code, x, y), dpi=300)


def plot_path(axes, ndf, x, y, vmin, vmax, cmap):
    l = ndf.shape[0]
    edf = pd.DataFrame({'i': np.arange(l-1), 'j': np.arange(1, l)})
    mplot.plot_nodes(axes, ndf, x, y, size=40, cmap=cmap, cbar=False,
                     zorder=4, lw=1.5, vmax=vmax, vmin=vmin)
    mplot.plot_edges(axes, ndf, edf, x, y, width=1.5, alpha=1, zorder=3,
                     color='black')


def read_landscape(code):
    fpath = 'output/{}.filtered.nodes.pq'.format(code)
    ndf = read_dataframe(fpath)
    fpath = 'output/protein.nodes.pq'
    protein = read_dataframe(fpath)[['region']]
    ndf = ndf.join(protein, on='protein')
    
    fpath = 'output/{}.filtered.edges.npz'.format(code)
    edf = read_edges(fpath)
    return(ndf, edf)


if __name__ == '__main__':
    vmin, vmax = (-8.078194290399551, 2.297213241457939)
    cmap = 'viridis'
    wt = 'VDGV'
    suffixes = ['FA', 'FG', 'FT', 'FC', 'FV',
                'LA', 'LG', 'LT', 'LC', 'LV',
                'CC', 'CA', 'CG',
                'SA', 'AA', 'AC', 'AG',
                'GA', 'GC', 'GV', 'GL', 'GM', 'GF', 'GT', 'GG', 'GV']

    code, x, y = 'Standard', '2', '1'
    code, x, y = sys.argv[1:]
    
    print('Reading filtered landscape under code {}'.format(code))
    ndf, edf = read_landscape(code)
    make_plot(ndf, edf, x, y, code, suffixes, size=5, padding=0.1)
