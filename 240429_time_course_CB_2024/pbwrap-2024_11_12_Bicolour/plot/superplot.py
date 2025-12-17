import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.container as container
import pbwrap.plot.violin as violin
import functools

from .utils import make_color_frame, get_markers_generator, get_colors_list

def distribution_super_plot(data, ax, title=None, xlabel=None, ylabel=None, showextrema=False) :
    """
    Distribution plot for pd.Series data.
    Series is expected to have a multi-index from 1 to 3 dimensions and Series values are expected to be lists.

    Super plot is returned as plt.Axes

    """
    if type(data) != pd.Series : raise TypeError('data argument passed is of type {0} should be "pd.Series".'.format(type(data)))
    if type(ax) != plt.Axes : raise TypeError('ax argument passed is of type {0} should be "plt.Axes".'.format(type(ax)))

    level = len(data.index.names)

    if level == 1 :
        ax, legend = _distribution_lvl1(data, ax, showextrema=showextrema)
    
    elif level == 2 :
        ax, legend = _distribution_lvl2(data, ax, showextrema=showextrema)

    elif level == 3 :
        ax, legend = _distribution_lvl3(data, ax, showextrema=showextrema)
    
    else : raise LevelError("Unsupported number of dimension in index (should be between 1 and 3).")
    
    legend_size = len(legend[0])
    if legend_size > 5 :
        ncols= (legend_size // 5) + 1
        legend = ax.legend(*legend, ncols=ncols)
    else : legend = ax.legend(*legend)
    
    if type(title) != type(None) : ax.set_title(title)
    if type(xlabel) != type(None) : ax.set_xlabel(xlabel)
    if type(ylabel) != type(None) : ax.set_ylabel(ylabel)

    return ax


class LevelError(IndexError) :
    pass

def _distribution_lvl1(data: pd.Series, ax: plt.Axes, showextrema=False) :

    multi_index = data.index.names
    if len(multi_index) != 1 : raise LevelError("_distribution_lvl1 was called but multi-index dimension does not match.")
    if type(data.iat[0]) != list : raise ValueError("for distribution plot a list is expected for each distributions to plot. {0} was found.".format(type(data.iat[0])))

    colors = get_colors_list(len(data))
    labels = list(data.index)

    ax = violin.violin_plot(
        ax=ax,
        distributions=data,
        labels= labels,
        colors= colors,
        showmeans= True,
        mean_size= 90,
        showextrema=showextrema,
    )

    legend_prep =(
        [plt.scatter(0,0, c= 'white', edgecolors= 'black')],
        ['distribution means']
        )

    return ax, legend_prep

def _distribution_lvl2(data: pd.Series, ax : plt.Axes, alpha= 0.6, showextrema=False, show_distribution_size=True) :
    
    multi_index = data.index.names
    if len(multi_index) != 2 : raise LevelError("_distribution_lvl2 was called but multi-index dimension does not match.")
    data = data.sort_index()
    measure = data.name
    distributions: pd.Series = data.reset_index(drop=False).groupby(multi_index[:1])[measure].apply(list)
    labels_lvl2 = list(distributions.index)
    if show_distribution_size :
        for index, label in enumerate(labels_lvl2) : 
            labels_lvl2[index] = label + '\n' + str([len(dis) for dis in distributions.iat[index]])
    labels_lvl1 = list(data.index.get_level_values(1).unique())

    colors_df = make_color_frame(labels= labels_lvl1)

    colors = list(pd.merge(data.reset_index(), colors_df, left_on= multi_index[1], right_on='labels').sort_values(multi_index)['colors'])

    ax = violin.violin_plot(
        ax=ax,
        distributions=distributions,
        labels= labels_lvl2,
        colors=colors,
        alpha=alpha,
        linewith=1,
        showmeans= True,
        mean_size= 90,
        showextrema=showextrema,
        multi_violin_plot= True
    )

    legend_prep = (
        [plt.scatter(0,0, c= 'white', edgecolors= 'black')] + [plt.bar(0,0, color= color) for color in colors_df['colors']],
        ['distribution means'] + list(colors_df.index)        
    )


    return ax, legend_prep


def _extract_color_from_legend(legend_prep, remove_first=False) :
    if remove_first :
        legend_prep = (legend_prep[0][1:], legend_prep[1][1:])
    
    artists: list[container.BarContainer] = legend_prep[0]
    colors = [(artist.patches[0].get_facecolor(),) for artist in artists]

    color_df = pd.DataFrame(data= colors, index= legend_prep[1], columns= ['colors'])
    return color_df

def _distribution_lvl3(data: pd.Series, ax: plt.Axes, showextrema) :
    """

    yaxis = measure_value

    level1 = x-axis (axis label)
    level2 = distribution colors (legend)
    level3 = smaller scatter points within distribution whith different shapes (legend)

    """
    multi_index = data.index.names
    if len(multi_index) != 3 : raise LevelError("_distribution_lvl3 was called but multi-index dimension does not match.")
    data = data.sort_index()
    measure = data.name

    #Levevel 1&2
    list_flattener = functools.partial(sum, start=[])
    data_distribution_lvl2: pd.Series = data.reset_index(drop=False).groupby(multi_index[:2])[measure].apply(list).apply(list_flattener)
    ax, legend_lvl2 = _distribution_lvl2(
        data_distribution_lvl2, 
        ax,
        alpha=0.4,
        showextrema=showextrema
        )

    #Level 3
    index_lvl3 = data.index.get_level_values(2).unique()

    gen = get_markers_generator()
    next(gen) #skipping 'o' that is used for means of distributions
    marker_list = [next(gen) for i in range(len(index_lvl3))]
    marker_frame = pd.Series(data= marker_list, index= index_lvl3, name= 'markers')
    color_df = _extract_color_from_legend(legend_lvl2, remove_first= True) #First artist is always distributions mean.

    markers = pd.merge(data.reset_index(), marker_frame, left_on= multi_index[2], right_on=multi_index[2]).sort_values(multi_index)['markers']
    colors = pd.merge(data.reset_index(), color_df, left_on= multi_index[1], right_index=True ).sort_values(multi_index)['colors']

    # positions, _ = violin.multi_violin_plot_positions(data_distribution_lvl2.reset_index().groupby(multi_index[0])[measure].apply(list))
    positions = [[i+1]*len(data.xs(idx,axis=0,level=1)) for i,idx in enumerate(data.index.get_level_values(1).unique())]
    print(positions)
    positions = sum(positions,[])
    # positions *= len(index_lvl3)
    # positions.sort()

    #random span for points
    seed = 0
    random_gen = np.random.default_rng(seed=seed)
    step_between_positions = 0.20

    assert len(markers) == len(colors) == len(positions) == len(data), 'all length should be identical :\nmarkers length : {0}\ncolors length : {1}\npositions length {2}\ndata length {3}'.format(len(markers), len(colors), len(positions), len(data))

    for position, distribution_lvl3, marker, color in zip(positions, data, markers, colors) :
        x = np.array([position] * len(distribution_lvl3), dtype=float)
        x += (random_gen.random(size=len(distribution_lvl3), dtype=float) * step_between_positions * random_gen.choice([1,-1], size=len(distribution_lvl3)))
        ax.scatter(
            x= x,
            y= distribution_lvl3,
            s = 25,
            color = color,
            marker = marker,
            alpha= 0.1,
            zorder =-100
        )

    legend_prep = (
        [plt.scatter(0,0, c= 'white', edgecolors= 'black')]  + [plt.bar(0,0, color= color) for color in color_df['colors']],
        ['distributions means'] + list(color_df.index)
    )

    return ax, legend_prep
