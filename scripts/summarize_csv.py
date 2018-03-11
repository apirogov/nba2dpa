#!/usr/bin/env python3
import pandas as pd
import seaborn as sns

import sys
df = pd.read_csv(sys.argv[1])

df.tool = df.tool.str.replace(r'cat .*ltl_','').str.replace('.sh > %O','')
tools = ['dstar','spot','rabinizer','nbadet']
clrs = {'dstar': 'r', 'spot': 'g', 'nbadet': 'b', 'rabinizer':'k'}

fs = list(set(df.formula))
fs = sorted(fs, key=lambda f: df.query("formula == '"+f+"' & tool=='nbadet'").states.sum())

ys = {}
for tool in tools:
    view = df[df.tool==tool]
    ys[tool] = []
    for f in fs:
        ys[tool].append(view[view.formula==f].states.sum())
    ax = sns.pointplot(list(range(len(fs))), ys[tool], color=clrs[tool], join=False)
    ax.set(yscale="log")
sns.plt.show()
# print(tool)
# print(df[df.tool == tool].states.describe())
