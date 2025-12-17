import matplotlib.pyplot as plt
import pandas as pd

def add_significance(ax: plt.Axes, pvalues: pd.DataFrame, significance = 0.01) :
    """
    Add a star in the label.
    """
    
    xtickslabels = ax.get_xticklabels()

    if len(xtickslabels) != len(pvalues) : raise ValueError("Expected same length for pvalues and x_ticks.\n pvalues : {0} x_ticks : {1}".format(len(pvalues), len(xtickslabels)))
    if isinstance(pvalues,pd.DataFrame) :
        if len(pvalues.columns) != 1 : raise ValueError("Pass Dataframe with loc on pvalues (1 col expected)")
    elif isinstance(pvalues, pd.Series) :
        pass
    else :
        raise TypeError("Expected dataframe or Series")
    
    res_labels = []

    for label, pvalue in zip(xtickslabels, pvalues) :
        label = label.get_text()
        if pvalue <= significance :
            res_labels.append(
                label + "\n*"
            )
        else :
            res_labels.append(label)
    
    ax.set_xticklabels(res_labels)

    return ax


def plot_pvalue_bar(ax: plt.Axes, x_sample1, y_sample1, x_sample2, y_sample2) :

    y_ticks = ax.get_yticks()
    if len(y_ticks) >= 2 : y_step = y_ticks[1] - y_ticks[0]
    else : y_step = 0

    y_sample1 +=y_step*0.2
    y_sample2 +=y_step*0.2

    ceiling = max(y_sample1 + y_step , y_sample2 + y_step,)

    ax.plot([x_sample1 - 0.1, x_sample1 + 0.1, x_sample1,x_sample1,x_sample2, x_sample2, x_sample2 - 0.1, x_sample2 + 0.1], [y_sample1, y_sample1, y_sample1, ceiling, ceiling, y_sample2, y_sample2, y_sample2], 'r')
    ax.text(
        x = x_sample1 + (x_sample2 - x_sample1)/2,
        y= ceiling,
        s= "*",
        c='red',
        fontdict={'size' : 20}
        )
    
    ax.set_ylim(top= ceiling)

    return ax


def add_pvalue_star(ax : plt.Axes) :
    """
    Add a star in the title.
    """
    title = ax.get_title()
    title += "\n*"
    ax.set_title(title)
    return ax