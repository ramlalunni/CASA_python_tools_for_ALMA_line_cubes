import numpy as np

# function to change axes colour
###----------------------------------------------------------------------------###
def ChangeAxisColor(ax,color):
    """Changes color of input axis and ticks (both major and minor) to input color"""

    sides = ['bottom','top','right','left']
    for side in sides:
        ax.spines[side].set_color(color)
        ax.tick_params(axis='x', colors=color,which='major',length=6,labelcolor="black", direction='in')
        ax.tick_params(axis='x', colors=color,which='minor',length=3,labelcolor="black", direction='in')
        ax.tick_params(axis='y', colors=color,which='major',length=6,labelcolor="black", direction='in')
        ax.tick_params(axis='y', colors=color,which='minor',length=3,labelcolor="black", direction='in')
    return ax
###----------------------------------------------------------------------------###


# function to align zeros of twin y axes
###----------------------------------------------------------------------------###
def align_yaxis_origins(axes):
    """Aligns origin point of the two y-axes created by matplotlib.axes.Axes.twinx()
       Taken from: https://stackoverflow.com/a/46901839/ (Licensed as CC BY-SA 4.0: https://creativecommons.org/licenses/by-sa/4.0/)
       See https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.twinx.html for Axes.twinx() documentation"""

    # converting input to numpy array
    axes = np.array(axes)

    # getting y axis limits
    extrema = np.array([ax.get_ylim() for ax in axes])

    # reset for divide by zero issues
    for i in range(len(extrema)):
        if np.isclose(extrema[i, 0], 0.0):
            extrema[i, 0] = -1
        if np.isclose(extrema[i, 1], 0.0):
            extrema[i, 1] = 1

    # upper and lower limits
    lowers = extrema[:, 0]
    uppers = extrema[:, 1]

    # if all pos or all neg, don't scale
    all_positive = False
    all_negative = False
    if lowers.min() > 0.0:
        all_positive = True

    if uppers.max() < 0.0:
        all_negative = True

    if all_negative or all_positive:
        # don't scale
        return

    # pick "most centered" axis
    res = abs(uppers+lowers)
    min_index = np.argmin(res)

    # scale positive or negative part
    multiplier1 = abs(uppers[min_index]/lowers[min_index])
    multiplier2 = abs(lowers[min_index]/uppers[min_index])

    for i in range(len(extrema)):
        # scale positive or negative part based on which induces valid
        if i != min_index:
            lower_change = extrema[i, 1] * -1*multiplier2
            upper_change = extrema[i, 0] * -1*multiplier1
            if upper_change < extrema[i, 1]:
                extrema[i, 0] = lower_change
            else:
                extrema[i, 1] = upper_change

        # bump by 10% for a margin
        extrema[i, 0] *= 1.1
        extrema[i, 1] *= 1.1

    # set axes limits
    [axes[i].set_ylim(*extrema[i]) for i in range(len(extrema))]
###----------------------------------------------------------------------------###
