#!/usr/bin/env python3
import pandas as pd
#import seaborn as sns

import sys
df = pd.read_csv(sys.argv[1])

df = df[['input.name','tool','exit_status','time','output.states']]
df.tool = df.tool.str.replace(r'cat .*\| ','').str.replace(r'.*/bin/','').str.replace(' > %O','')

tools = df.tool.unique()
samples = df.groupby('input.name')
inputs = samples.groups

if (len(sys.argv)>2):
    inputs = inputs[0:int(sys.argv[2])]

print(len(tools),"tools,",len(inputs),"inputs")

#print("determining complete samples...")
goodinputs = []
for name, grp in samples:
    if len(grp[grp.exit_status=='ok']) == len(tools):
        goodinputs.append(name)
goodinputs = set(goodinputs)

#print("removing invalid samples...")
df = df[df['input.name'].isin(goodinputs)]

print(len(goodinputs),"complete samples. total states per tool:")
toollist = df.groupby('tool')
grps = sorted(toollist.groups)
refsum = toollist.get_group(grps[0])['output.states'].sum()
for t in grps:
    statesum = toollist.get_group(t)['output.states'].sum()
    print(t, statesum, "{0:.2f}".format(statesum/refsum))
