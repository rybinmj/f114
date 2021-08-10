import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit as cf
import numpy as np
from statannot import add_stat_annotation as anot

class Curve():
    def __init__(self):
        self.rc = {
            'ytick.labelsize': 16,
            'xtick.labelsize': 18,
            'font.size': 25,
        }
        self.colors = sns.color_palette('muted').as_hex()

        
def four_parameter_logistic(x, A, B, C, D):
    return ((A-D)/(1.0+((x/C)**(B))) + D)


def run_4pl(df):
    m = df.melt()
    x = [float(val) for val in m['variable'].values.tolist()[6:]]
    y = m['value'].tolist()[6:]
    params, params_covariance = cf(four_parameter_logistic, x, y)
    return params, params_covariance
    

def plot_legend(params, color='black', label='default_label'):
    ic50 = u'IC\u2085\u2080'
    ic50_value = str("{:.2f}".format(params[2]*1000000))
    patch = mpatches.Patch(color=color, label=(f'{label}, {ic50}: {ic50_value} µM'))
    return patch


def plot_lines(df, color='black'):
    m = df.melt()
    x = [float(val) for val in m['variable'].values.tolist()[6:]]
    y = m['value'].tolist()[6:]

    params, cov = run_4pl(df)
    x_line = np.arange(0.0, max(x)*1.1, min(x))
    y_line = four_parameter_logistic(x_line, *params)
    ax = sns.lineplot(
        x=x_line, y=y_line, linewidth=2, solid_joinstyle='round', color=color,
        err_style='band', 
        zorder=20,
    )
    ax = sns.lineplot(
        x=x, y=y, 
        color=color, marker='.', markersize=20, linestyle='', 
        err_style='bars', ci=68, err_kws={'capsize':7, 'elinewidth':2, 'capthick':2}, 
        zorder=0,
    )
    return params

    
def fit(
    data,
    rc=Curve().rc, 
    figsize=[6.5,6],
    ylim=[-10, 155],
    legend_loc='lower left',
    ylabel='% Control',
    xlabel='log([F114]) M',
    **kwargs
):
    if type(data) != list:
        data = [data]
        
    plt.figure()
    plt.rcParams.update(**rc)
    fig, ax = plt.subplots(figsize=[6.5,6])

        
    patches = []
    for n, df in enumerate(data):
        c = kwargs['color'][n]
        l = kwargs['label'][n]
        params = plot_lines(df, color=c)
        patches.append(plot_legend(params, color=c, label=l))

    sns.despine()
    ax.set_xscale('log')
    ax.set_ylim(ylim)
    plt.legend(handles=patches, frameon=False, loc=legend_loc, fontsize=15)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel, labelpad=6)
    plt.tight_layout()
    
    

class Bar:
    def __init__(self):
        self.rc = {
            'font.size': 20.0,
            'axes.labelsize': 20.0,
            'axes.titlesize': 'large',
            'xtick.labelsize': 20,
            'ytick.labelsize': 16,
            'legend.fontsize': 'medium',
            'font.family': 'Arial',
        }
        
        self.colorblind_palette = sns.color_palette('colorblind').as_hex()
        self.dark_palette = sns.color_palette('dark').as_hex()

        self.grey_and_green = [
            self.colorblind_palette[7],
            self.dark_palette[2],
        ]
        
def plot_bargraph(
    df, 
    color=Bar().grey_and_green, 
    xticklabels=['DMSO', 'F114'], ylabel='% GFP+ Cells (Normalized to DMSO)', 
    show_thresholds=1
):
    plt.figure(figsize=[6,6])
    plt.rcParams.update(**Bar().rc)
    plot = sns.barplot(data=df, palette=color, errcolor='black', capsize=0.1, errwidth=3, ci=68)
    sns.stripplot(data=df, color='black', size=8)
    plt.xticks(ticks=[0,1], labels=xticklabels)
    plt.ylabel(ylabel)
    anot(
        plot, data=df,
        test='t-test_ind', box_pairs=[('dmso', 'f114')], comparisons_correction=None,
        fontsize=24, line_offset_to_box=0.1, line_height=0.05, text_offset=-6, color='black', verbose=show_thresholds,
    )
    sns.despine()
    plt.tight_layout()