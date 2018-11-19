#!/usr/bin/env python3
import pandas as pd
#import seaborn as sns

import sys
df = pd.read_csv(sys.argv[1])

df = df[['input.name','tool','exit_status','time','output.states']]
df.tool = df.tool.str.replace(r'cat .*\| ','').str.replace(r'.*/bin/','').str.replace(' > %O','')

tools = df.tool.unique()
inputs = df['input.name'].unique()

if (len(sys.argv)>2):
    inputs = inputs[0:int(sys.argv[2])]

print(len(tools),"tools,",len(inputs),"inputs")

def numok(df, input):
    return len(df[df['input.name']==input][df.exit_status=='ok'])

print("determining complete samples...")
goodinputs = []
for inp in inputs:
    if numok(df, inp) == len(tools):
        goodinputs.append(inp)
goodinputs = set(goodinputs)

df = df[df['input.name'].isin(goodinputs)]

print("calculate total sums...")
sts = []
for tool in tools:
    sts.append(df[df.tool == tool]['output.states'].sum())

print(len(goodinputs),"complete samples. total states per tool:")
for i in range(len(tools)):
    print(tools[i], sts[i])
