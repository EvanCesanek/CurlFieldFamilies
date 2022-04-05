from cmath import nan
import matplotlib.pyplot as plt
import pandas as pd  # for data analysis
from scipy import stats
import numpy as np  # for working with arrays
import json  # for loading JSON format data
import os  # for navigating the file system
import seaborn as sns
custom_params = {'axes.spines.right': False,
                 'axes.spines.top': False,
                 'figure.dpi': 300,
                 'savefig.dpi': 300}
sns.set_theme(context='talk',  # paper, notebook, talk, poster
              style='ticks',
              palette=sns.color_palette('deep', 20),
              rc=custom_params)
object_colormap = 'OrRd'

x,y = np.meshgrid(np.linspace(-5,5,10),np.linspace(-5,5,10))

fig, ax = plt.subplots(1,5,figsize=(10,2.6), sharey=True, dpi=300)
for i in np.arange(0,5):
    scale = .015+0.005*i
    if i==2: scale = .015+0.005*4
    u = y*scale
    v = -x*scale
    ax[i].quiver(x,y,u,v,scale=1)
ax[2].set_xlabel('X Velocity')
ax[0].set_ylabel('Y Velocity')
fig.tight_layout()
fig.savefig('curl_fields.pdf')