"""
Collects together convenient plotting routines.
"""
import matplotlib.pyplot as plt
import seaborn as sns
import math
import numpy


def mm_to_in(mm):
    return float(mm) / 25.4


def plot_histogram(data, filename, xlabel, ylabel, xLim=None, title=None):
    """ Generic method for plotting a histogram. """
    if not data:
        # list is empty!
        print('WARNING: cannot create plot ' + filename + ' as no data was supplied to plot.')
        return
    data = [d for d in data if not math.isnan(d)]   # remove 'nan' from data.
    plt.clf()
    plt.hist(data)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if xLim:
        plt.xlim(xLim)
    if title:
        plt.title(title)
    font = {'size': 18}
    plt.rc('font', **font)
    plt.savefig(filename + '.png', dpi=300)
    plt.savefig(filename + '.eps', dpi=300)


def ecdf(observations):
    """
    Empirical Cumulative Distribution Function.
    :param observations: list of observations, in any order.
    :return: a list of tuples: [ (<value>,<proportion of observations less than or equal to value>) ]
    Values will be given in ascending order.
    """
    s = sorted(observations)  # non-destructive, preserves original order of observations.
    n = len(s)  # number of observations.
    cdf = []
    for i, sample in enumerate(s):
        proportion = float(i + 1) / n
        tup = (sample, proportion)
        cdf.append(tup)
    return cdf


def plot_cdfs(
        distributions: list,
        vertical_xs=[],  # List of values at which to vertical lines indicating values of note.
        names=None,
        colours=None,
        filename=None,
        # If false, use seaborn to plot the distribution
        empirical=True,
        hist=True, rug=True, cumulative=False,
        # Other options
        show_legend=True,
        xlabel=None, ylabel=None, title=None, despine_y=False,
        fontsize=7,
        linewidth=1.0,
        # X axis manipulation
        xlog=False,  # False, linear. True, log scale. 'symlog', symetrical log scale that has a section of linear around zero.
        x_symlog_thresh=None,  # Float, if used.
        xmin=None, xmax=None,
        line_styles=None,
        plot_uniform=False,
        figsize_mm=(100, 100)  # In mm
        ):
    """
    Calculates cumulative distribution functions for each distribution in the supplied list.
    :param distributions: list of lists. Each lower level list is a distribution to plot as empirical cdf
    :param names: name of each distribution, used for legend. List of strings.
    :param colours: list of objects interpretable as colours. If used, one per distribution. E.g. ['r','k','b','g']
    :param xlabel: xlabel of graph
    :param title: title of graph
    :param filename: if left as None, will show graph visually rather than save to filesystem. Leave off the extension
    (.png and similar), as these are added. Files saved as both eps and png.
    :param fontsize:
    :param linewidth:
    :param xlog: False, linear scale. True, log scale. 'symlog', symetrical log scale that has a section of linear
    around zero.
    :param x_symlog_thresh: Thresholds around zero for linear plot if 'symlog' selected.
    :param plot_uniform: if True, plot the distribution for a uniform distribution within the max and min values of the
    other distributions shown.
    """
    f, ax = plt.subplots(figsize=[mm_to_in(figsize_mm[0]), mm_to_in(figsize_mm[1])],
                         constrained_layout=True)
    font = {'size': fontsize}
    plt.rc('font', **font)
    # sns.set_style("darkgrid")
    plots = []
    if line_styles is None:
        line_styles = ['-'] * len(distributions)
    legend_names = []
    if names is None:
        names = [None] * len(distributions)
    for i, (distro, name, ls) in enumerate(zip(distributions, names, line_styles)):
        if len(distro) > 0:
            col = colours[i] if colours else None
            if empirical:
                cdf = ecdf(distro)
                # insert one record at start of these lists to ensure curves extend to y=0.
                D1X = [cdf[0][0]]
                D1X.extend([record[0] for record in cdf])
                D1Y = [0]
                D1Y.extend([record[1] for record in cdf])

                p1, = ax.plot(D1X, D1Y, linewidth=linewidth, color=col)
                plots.append(p1)
                legend_names.append(name)
            else:  # KDE smoothed. Does not work for distributions containing 1 item, and can be misleading.
                sns.distplot(distro, ax=ax, hist=hist, kde=True, rug=rug, norm_hist=True,
                             kde_kws=dict(cumulative=cumulative), label=name, color=col)

    # Vertical line indicating points of note
    for vi, vx in enumerate(vertical_xs):
        if colours:
            plt.plot((vx, vx), (0., 1.), linestyle=':', color=colours[vi])
        else:
            plt.plot((vx, vx), (0., 1.), linestyle=':')

    if plot_uniform:
        all_data = numpy.concatenate([numpy.array(d) for d in distributions])
        x_min = numpy.min(all_data)
        x_max = numpy.max(all_data)
        plt.plot((x_min, x_max), (0., 1.), linewidth=linewidth)
    if len(legend_names) > 0 and None not in legend_names and show_legend:
        plt.legend(plots, legend_names, loc='lower right')
    plt.ylabel(ylabel)
    if xlabel:
        plt.xlabel(xlabel)
    if title:
        plt.title(title)
    plt.gca().grid(True, linewidth=linewidth)   # turn on grid lines.
    # Change the width of the plot boundary lines
    ax = plt.gca()
    [i.set_linewidth(linewidth) for i in ax.spines.values()]
    if xmax is not None or xmin is not None:
        curXmin, curXmax = plt.gca().get_xlim()
        if xmax is not None:
            curXmax = xmax
        if xmin is not None:
            curXmin = xmin
        plt.gca().set_xlim([curXmin, curXmax])
    if xlog is True:
        plt.xscale('log')
    if xlog is 'symlog':
        if x_symlog_thresh:
            plt.xscale('symlog', linthreshx=x_symlog_thresh)
        else:
            plt.xscale('symlog')
    if cumulative or empirical:
        pass
    else:
        sns.despine(left=True)
    if despine_y:
        plt.tick_params(
            axis='y',  # changes apply to the y-axis
            which='both',  # both major and minor ticks are affected
            left='off',  # ticks along the bottom edge are off
            right='off',  # ticks along the top edge are off
            labelleft='off')  # labels along the left edge are off
    if filename:
        plt.savefig(filename + '.pdf')  # , transparent=True)  # Transparent doesn't work well if using sns dark grid
        plt.savefig(filename + '.png', dpi=600)
        # plt.savefig(filename + '.eps', dpi=600)
    else:
        plt.show()


# Testing
# dist = [1,2,2,234,6,7,3,23,214,7,6,83,312,32,65,78,58]
# plot_cdfs([dist], figsize_mm=(50, 50), filename='test', empirical=True)