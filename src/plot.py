import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# Define Colors
SEEBLAU = '#00A9E0'
SEEGRAU = '#9AA0A7'
SEEBLAU_LIST = ['#CCEEF9', '#A6E1F4', '#59C7EB', '#00A9E0', '#008ECE']


fontsize = 10
fontsize_small = 9
figsize = np.array((6, 2.5))
figsize_small = np.array((3, 1.8))

# reset matplotlib
mpl.rcdefaults()

# change dpi
plt.rcParams.update({
	'savefig.dpi': 300,
	'figure.dpi': 150,
	'figure.figsize': figsize,
	'figure.autolayout': True,
    'savefig.bbox': 'tight',
})

# change font type
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': 'Roboto',
    'text.usetex': False,
})

# change font sizes
plt.rcParams.update({
    'font.size': fontsize,
    'axes.titlesize': fontsize,
    'axes.labelsize': fontsize,
    'xtick.labelsize': fontsize_small,
    'ytick.labelsize': fontsize_small,
    'legend.fontsize': fontsize,
    'figure.titlesize': fontsize_small,
    'axes.labelpad': 2,
    'axes.labelsize': fontsize_small,
})

plt.rcParams.update({
    'errorbar.capsize': 2,
})

# setup ticks
plt.rcParams.update({
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'xtick.minor.size': 2.5,
    'ytick.minor.size': 2.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.width': 0.5,
    'ytick.minor.width': 0.5,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.top': True,
    'ytick.right': True,
})

# change the color of the box and the ticks
plt.rcParams.update({
    'axes.edgecolor': SEEGRAU,
    'xtick.color': SEEGRAU,
    'ytick.color': SEEGRAU,
    'xtick.labelcolor': 'black',
    'ytick.labelcolor': 'black',
})

# setup grid
plt.rcParams.update({
    'axes.grid': True,
    'grid.color': SEEGRAU,
    'grid.linewidth': 0.5,
    'grid.alpha': 0.2,
})

# setup cmaps
cmap = 'magma'
cmap_diverging = 'seismic'
cmap_cyclic = 'twilight'
plt.rcParams['image.cmap'] = cmap

# setup colors
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=
[(0.0, 0.6627450980392157, 0.8784313725490196),
 (0.6039215686274509, 0.6274509803921569, 0.6549019607843137),
 (0.24313725490196078, 0.32941176470588235, 0.5882352941176471),
 (0.5568627450980392, 0.12549019607843137, 0.2627450980392157),
 (0.8784313725490196, 0.3764705882352941, 0.49411764705882355),
 (0.996078431372549, 0.6274509803921569, 0.5647058823529412),
 (0.0392156862745098, 0.5647058823529412, 0.5254901960784314),
 (0.027450980392156862, 0.44313725490196076, 0.5294117647058824)]
)
