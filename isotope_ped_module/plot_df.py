# -*-Python-*-
# Created by chabanr at 14 Dec 2021  13:49

"""
This script <fill in purpose>

defaultVars parameters
----------------------
:param kw1: kw1 can be passed to this script as <path to script>.run(kw1='hello')
"""
import matplotlib as mpl

testc = ['PTOT', 'PTOT']
testy = ['PRMPED_NESEP', 'PRMPED_TESEP']
testx = ['PRMPED_NEPED', 'PRMPED_TEPED']

# defaultVars(abscissa=['PRMPED_NEPED'], ordinate=['PRMPED_NESEP'], colorby=['PTOT'])
defaultVars(abscissa=testx, ordinate=testy, colorby=testc, logcolor=False, errors=False, colormap='jet', df=None)
for Input in [abscissa, ordinate, colorby]:
    if type(Input) not in [list, tuple]:
        print("All inputs must be iterables")

if df is None:
    df = root['OUTPUTS']['df']
# df_err = root['OUTPUTS']['df_err']

# check that they are all the same length
if all(len(item) == len(abscissa) for item in [abscissa, ordinate, colorby]):
    print("all input arrays same length all clear to continue")
# elif len(colorby) == 1:  # TODO: Look at making them all the same color using a single colorby


# f, ax = plt.subplots(nrows=len(abscissa) // 4 + 1, ncols=min([len(abscissa), 4]))
# try:
#     ax = ax.flatten()  # if it cant be flattened then its the single axes object
# except AttributeError:
#     ax = [ax]

# Create new figure and axis for each plot
f, ax = [], []
for i in range(len(abscissa)):
    fig, axes = plt.subplots()
    ax.append(axes)
    f.append(f)


# separate the isotopes:
dfH = df.where(df['ISOTOPE'] == 'H')
dfD = df.where(df['ISOTOPE'] == 'D')

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

i = 0
for xx, yy, cc in zip(abscissa, ordinate, colorby):

    # create normalization
    cmin, cmax = np.nanmin(df[cc].values), np.nanmax(df[cc].values)
    if logcolor:  # use a logarithmic colorbar
        cnorm = mpl.colors.LogNorm(vmin=cmin, vmax=cmax)
    else:
        cnorm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
    # print(np.min(df[cc].values), np.max(df[cc].values))
    cmap = mpl.cm.get_cmap(colormap)

    ms = ['s', 'o']  # marker style
    for q, dfX in enumerate([dfH, dfD]):
        dfX = dfX.dropna(subset=[xx, yy, cc])  # drop rows with missing values
        # print(dfX.head())

        cdat = cnorm(dfX[cc].values)  # sample the colormap
        # print(color_samples)
        dfX.plot.scatter(xx, yy, c=cmap(cdat), marker=ms[q], s=75, ax=ax[i], zorder=2)  # plot the markers

        if errors:
            ax[i].errorbar(
                dfX[xx], dfX[yy], xerr=dfX[xx + '_ERR'], yerr=dfX[yy + '_ERR'], fmt='o', marker=None, ms=0, mew=0, c='k', zorder=1
            )  # plot the errorbars
    # ENDFOR dfX

    # cbaxes = inset_axes(ax[i], width="5%", height="80%", loc=1) # used when I want to put colorbar inside axis
    cm_mappable = mpl.cm.ScalarMappable(norm=cnorm, cmap=colormap)
    # plt.colorbar(mappable=cm_mappable, cax=cbaxes, ticks=np.linspace(cmin, cmax, 5), orientation='vertical') # places colorbar inside axis
    plt.colorbar(
        mappable=cm_mappable, ax=ax[i], orientation='vertical'  # ticks=np.linspace(cmin, cmax, 5), orientation='vertical'
    )  # steals room for colorbar from axis
    ax[i].set_title(f"Color: {cc}")
    i += 1
