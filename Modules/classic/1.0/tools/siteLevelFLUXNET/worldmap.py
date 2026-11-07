import pandas as pd
from pandas import DataFrame
import yaml
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import colorsys
import seaborn as sns

import sys
import os

# first argument is the Processed_Sites directory
inputFiles=sys.argv[1]

# optional second argument is the output path
if len(sys.argv) > 2:
    outputFiles = sys.argv[2] + '/'
else:
    outputFiles = ''

def yaml_parse(file, att):
    with open(file, "r") as f:
      d = yaml.safe_load(f)
    return(d[att])



# import my list of sites.

#df
sites = []
for site in os.listdir(inputFiles):
    if os.path.isdir(inputFiles + "/" + site):
        sites.append(site)

data = []
for site in sites:
    print("starting " + site)
    path = inputFiles + "/" + site + "/siteinfo.yaml"
    row = []
    row.append(site)
    row.append(yaml_parse(path, "lat"))
    row.append(yaml_parse(path, "lon"))
    row.append(yaml_parse(path, "biome"))
    data.append(row)
df = DataFrame(data, columns=["Site","latitude","longitude","biome"])
print(df)
# Find the unique biomes we have
biomes=df.biome.unique()

# Find unique colours for the number of biomes
N = len(biomes)
palette = sns.color_palette(None, N)

# %%
# now plot
fig = plt.figure(figsize=[10, 8])
ax2 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
fig.subplots_adjust(bottom=0.1, top=0.9,
                    left=0.1, right=0.9)
ax2.coastlines(resolution='50m')
# Sites by biome
for n,b in enumerate(biomes):
    bio=df[df['biome'].str.match(b)]
    bio.plot.scatter(x='longitude', y='latitude',s=45,ax=ax2,label=b,color=palette[n],
             linewidth=1,  marker='o',alpha=0.9,edgecolors='black',
             transform=ccrs.PlateCarree()
             )

ax2.set_extent([-180,180,-60,90],crs=ccrs.PlateCarree())

#texts = [plt.text(df.longitude[i], df.latitude[i], df.Site[i], ha='center',
#            va='center',transform=ccrs.PlateCarree()) for i in range(len(df.Site))]

#ax2.set_boundary(circle, transform=ax2.transAxes)

#Improve the legend
handles, labels = ax2.get_legend_handles_labels()
legend = ax2.legend(handles[:], labels[:], title="Biome", prop = {'size':'large'},
          handletextpad=1, columnspacing=1,
          loc="lower left", ncol=1, frameon=True, fontsize=16)

for item in ([ax2.xaxis.label, ax2.yaxis.label] +
             ax2.get_xticklabels() + ax2.get_yticklabels()):
    item.set_fontsize(16)

plt.setp(legend.get_title(),fontsize='16')

# %%
plt.savefig(str(outputFiles + '/mapofFLUXNETsites.png'),dpi=450, bbox_inches='tight')
plt.show()
