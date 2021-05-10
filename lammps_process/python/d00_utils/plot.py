


def include_0_y_axis(fig, ax):
    if ax.get_ylim()[0]*ax.get_ylim()[1] <= 0:
        pass
    else:
        if ax.get_ylim()[0] > 0:
            ax.set_ylim(bottom=0)
        else:
            ax.set_ylim(top=0)
    return (fig, ax)
