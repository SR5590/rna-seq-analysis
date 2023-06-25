import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly
import sys
from Bio import SeqIO
import gzip
from statistics import mean
types=["A","B","C"]

def GeneCount(featureTable, ID):
    table=featureTable
    df1 =pd.read_csv(table, sep="\t", skiprows=1)
    df2=pd.read_csv("/mnt/g27prist/CMTD/annotation/hs/Gencode/GRCh37.p13/gencode.v19.chr_patch_hapl_scaff_noHaplNorPatch.annotation.gtf",
                sep="\t", skiprows=5, header=None)
    df3=pd.DataFrame()
    df3["Geneid"]=df2[8].str.split(";").str[0].str.lstrip(" gene_id ").str.rstrip('"').str.lstrip('"')
    df3["gene type"]=df2[8].str.split(";").str[2].str.lstrip(" gene_type ").str.rstrip('"').str.lstrip('"')
    df3=df3.drop_duplicates()
    df4=pd.merge(df3, df1,on="Geneid")
    new_columns=[]
    for a in list(df4):
        a=a.split("/")
        new_columns.append(a[len(a)-1])
    df4.columns= new_columns

    read_count=[]
    for a in list(df4)[7:]:
        read_count.append(len(df4[df4[a]>10]))
    pandas_dict={}
    pandas_dict["ID"]=ID
    pandas_dict["0.1"]=read_count[0]
    pandas_dict["0.2"]=read_count[1]
    pandas_dict["0.3"]=read_count[2]
    pandas_dict["0.4"]=read_count[3]
    pandas_dict["0.5"]=read_count[4]
    pandas_dict["0.6"]=read_count[5]
    pandas_dict["0.7"]=read_count[6]
    pandas_dict["0.8"]=read_count[7]
    pandas_dict["0.9"]=read_count[8]
    pandas_dict["1"]=read_count[9]
    return pandas_dict

df=pd.DataFrame(columns=["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1", "ID"])
for a in types:
    df=df.append(GeneCount("../05_featureCount/sub_count/"+a+"1" +"/FeatureCountTable", a+"1"), ignore_index=True)
    df=df.append(GeneCount("../05_featureCount/sub_count/"+a+"1" +"/FeatureCountTable", a+"2"), ignore_index=True)

data=[]
for a in types:
    xx=list(df[df["ID"].str.startswith(a)].mean())
    data.append(go.Scatter(x=["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"], y=xx, name=a))

fig = go.Figure(data=data)
fig.show()
