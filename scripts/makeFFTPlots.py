#!/usr/local/bin/python
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import seaborn as sns
import numpy as np
import pandas as pd
sp = "100"
lp = "100"
fname = "soft_pf0.2_sp" + sp + "_lp" + lp + "_condensed.density"
df=pd.read_csv(fname,delim_whitespace=True,header=None)
fig,ax = plt.subplots(1,2,figsize=(8,3))
#sns.heatmap(df,cmap=cm.viridis,ax=ax[0])
data = df.replace(0,1e-10)
data=data/data.sum().sum()
min_data = data.min().min()
if (min_data == 0):
    min_data = 1
max_data = data.max().max()
log_norm = LogNorm(vmin=min_data, vmax=max_data)
cbar_ticks = [10**i for i in 
              range(int(np.floor(np.log10(min_data))),
                    1 + int(np.ceil(np.log10(max_data))))]
sns.heatmap(data, norm=log_norm,
            cmap=cm.viridis, ax=ax[0],
            cbar_kws={"ticks": cbar_ticks})
fft_data=np.fft.fftshift(np.fft.fft2(df))
data = np.abs(fft_data)
# data=data/data.sum().sum()
min_data = data.min().min()
if (min_data == 0):
    min_data = 1
max_data = data.max().max()
log_norm = LogNorm(vmin=min_data, vmax=max_data)
cbar_ticks = [10**i for i in 
              range(int(np.floor(np.log10(min_data))),
                    1 + int(np.ceil(np.log10(max_data))))]
sns.heatmap(data, norm=log_norm,
            cmap=cm.viridis, ax=ax[1],
            cbar_kws={"ticks": cbar_ticks})
savename = "sp"+sp+"_lp"+lp
fig.savefig(savename+".png",dpi=300)
f=open(savename+"_fft_max.txt",'w')
f.write(str(np.max(data[data.shape[0]//2])))
f.close()
