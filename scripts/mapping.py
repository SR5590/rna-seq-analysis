import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly
import sys
from Bio import SeqIO
import gzip

ID={}
group=[]
total=[]
uniquemapped=[]
nonuniquemapped=[]

#Annotation_file="annotation.csv"
#Annotation_file=sys.argv[1]
Analyses_Folder="../"
#Analyses_Folder="./"
#with open(Annotation_file) as t:
#    for line in t:
#        if not line.startswith("ID"):
#            line=line.rstrip()
#            line=line.split("\t")
#            ID[line[0]]=line[2]
#            group.append(line[2])

#Annotation_file=pd.read_csv(sys.argv[1], sep="\t")
Annotation_file=pd.read_csv("sample_working.tsv", sep="\t")

df=pd.DataFrame()
for a in list(Annotation_file["working"]):
    with open(Analyses_Folder+"03_mapped/"+a+"/"+a+"Log.final.out") as s:
        stats={}
        stats["ID"]=a
        stats["type"]=Annotation_file[Annotation_file["working"]==a]["Type"].values[0]
        for line in s:
            line=line.rstrip()
            line=line.split("\t")
            if "Uniquely mapped reads number" in line[0]:
                stats["Unique mapped"]=int(line[1])
            elif "Number of input reads" in line[0]:
                stats["input reads"]=int(line[1])
            elif "Number of reads mapped to multiple loci" in line[0]:
                stats["multiple loci"]=int(line[1])


    with open(Analyses_Folder+"04_Bedtools_count/"+ a +"/"+a +"_exonCount.txt") as s:
        for line in s:
            line=line.rstrip()
            stats["exon_count"]=int(line)
    with open(Analyses_Folder+"04_Bedtools_count/"+ a +"/"+a +"_intergenicCount.txt") as s:
        for line in s:
            line=line.rstrip()
            stats["intergenic_count"]=int(line)
    with open(Analyses_Folder+"04_Bedtools_count/"+ a +"/"+a +"_intronCount.txt") as s:
        for line in s:
            line=line.rstrip()
            stats["intron_count"]=int(line)
    df=df.append(stats, ignore_index=True)



fig=go.Figure()
fig.update_xaxes(range=[0, df["input reads"].max() + df["input reads"].max()/20],nticks=10)
fig.add_trace(
    go.Bar(name='Total reads', y=df["ID"], x=df["input reads"],
           orientation="h",marker=dict(color="crimson")))
fig.add_trace(
    go.Bar(name='Unique mapped reads', y=df["ID"], x=df["Unique mapped"],
            orientation="h", marker=dict(color="darkblue")))
fig.add_trace(
    go.Bar(name='Exon mapped', x=df["exon_count"], y=df["ID"],
            orientation="h",marker=dict(color="darkgreen")))
fig.add_trace(
    go.Bar(name='Intron mapped', x=df["intron_count"], y=df["ID"],
            orientation="h",marker=dict(color="goldenrod")))
fig.add_trace(
    go.Bar(name='Intergenic mapped', x=df["intergenic_count"], y=df["ID"],
            orientation="h",marker=dict(color="dodgerblue")))

fig.update_layout(barmode='group', title=go.layout.Title(text="<b>Mapping Overview</b>",xref="paper",
                                                         x=0.5, y=0.95, font={"size":20}))

#list(range(1,2))
shapes=[]
y_pos=-0.5
for type in set(Annotation_file["Type"]):
    if y_pos==-0.5:
        y_pos=y_pos+len(Annotation_file[Annotation_file["Type"]==type]["ID"])
    elif y_pos>-0.5:
        shapes.append(go.layout.Shape(type="line",x0=0,y0=y_pos, x1=df["input reads"].max() + df["input reads"].max()/20, y1=y_pos,line=dict(dash="dot",color="black",width=2)))
        y_pos=y_pos+len(Annotation_file[Annotation_file["Type"]==type]["ID"])

fig.update_layout(shapes=shapes,xaxis=go.layout.XAxis(title=go.layout.xaxis.Title(text="Number of Reads",
                                                  font=dict(family="Courier New, monospace",size=18,color="#7f7f7f"))))


first_one=plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')

print("first plot done")


##########################################################################################
group_dict={}
reads_dict={}
print("counting reads")
for a in list(Annotation_file["working"]):
    group_dict[a]=Annotation_file[Annotation_file["working"]==a]["Type"].values[0]
    with gzip.open(Analyses_Folder+"02_trimmed/"+a+"/"+a+"_trimmed_1P.fastq.gz","rt") as s:
        read_number=0
        for read in SeqIO.parse(s, "fastq"):
            read_number=read_number+1
        reads_dict[a]=read_number

print("finished counting")
df1 =pd.read_csv(Analyses_Folder+"05_featureCount/FeatureCountTable", sep="\t", skiprows=1)
#df1 =pd.read_csv("test1", sep="\t", skiprows=1)
df2=pd.read_csv("/mnt/g27prist/CMTD/annotation/hs/Gencode/GRCh37.p13/gencode.v19.chr_patch_hapl_scaff_noHaplNorPatch.annotation.gtf",
                sep="\t", skiprows=5, header=None)
df3=pd.DataFrame()
df3["Geneid"]=df2[8].str.split(";").str[0].str.lstrip(" gene_id ").str.rstrip('"').str.lstrip('"')
df3["gene type"]=df2[8].str.split(";").str[2].str.lstrip(" gene_type ").str.rstrip('"').str.lstrip('"')
df3=df3.drop_duplicates()
df4=pd.merge(df3, df1,on="Geneid")
df4.columns = df4.columns.str.replace("../03_mapped/","")
df4.columns = df4.columns.str.replace("Aligned.sortedByCoord.out.mkdup.bam","")
df4.columns = df4.columns.str.split("/").str[0]




print(list(df4))

for key,value in reads_dict.items():
    df4["RPKM_"+key]=df4[key]/((value/1000000)*(df4["Length"]/1000))


#################################################################################################################
number_of_groups=len(group)
fig=go.Figure()
data=[]
max_reads=df4[list(df4)[7:7+(len(reads_dict))]].sum().max()
fig.update_xaxes(range=[0, max_reads + max_reads/20],nticks=10)
for a in set(df4["gene type"]):
    df3_sub=df4[df4["gene type"]==a]
    df3_sub = df3_sub[list(df4)[7:7+(len(reads_dict))]]
    name=str(a)
    x=list(df4)[7:7+(len(reads_dict))]
    y=list(df3_sub.sum())
    fig.add_trace(
        go.Bar(name=name, y=x, x=y, orientation="h"))

shapes=[]
y_pos=-0.5
for type in set(Annotation_file["Type"]):
    if y_pos==-0.5:
        y_pos=y_pos+len(Annotation_file[Annotation_file["Type"]==type]["ID"])
    elif y_pos>-0.5:
        shapes.append(go.layout.Shape(type="line",x0=0,y0=y_pos, x1=max_reads + max_reads/20, y1=y_pos,line=dict(dash="dot",color="black",width=2)))
        y_pos=y_pos+len(Annotation_file[Annotation_file["Type"]==type]["ID"])

fig.update_layout(shapes=shapes,barmode='stack',title=go.layout.Title(text="<b>Raw featureCount distribution</b>",xref="paper", x=0.5, y=0.95, font={"size":20}))

fig.update_layout(xaxis=go.layout.XAxis(title=go.layout.xaxis.Title(text="Raw feature counts",
                                                  font=dict(family="Courier New, monospace",size=18,color="#7f7f7f"))))

Second_one=plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')

print("second plot done")
###############################################################################################################

number_of_groups=len(group)
fig=go.Figure()
data=[]
max_reads=df4[list(df4)[7+(len(reads_dict)):]].sum().max()
fig.update_xaxes(range=[0, max_reads + max_reads/20],nticks=10)
for a in set(df4["gene type"]):
    df3_sub=df4[df4["gene type"]==a]
    df3_sub = df3_sub[list(df4)[7+(len(reads_dict)):]]
    name=str(a)
    x=list(df4)[7+(len(reads_dict)):]
    y=list(df3_sub.sum())
    fig.add_trace(
        go.Bar(name=name, y=x, x=y, orientation="h"))


y_pos=-0.5
shapes=[]
for type in set(Annotation_file["Type"]):
    if y_pos==-0.5:
        y_pos=y_pos+len(Annotation_file[Annotation_file["Type"]==type]["ID"])
    elif y_pos>-0.5:
        shapes.append(go.layout.Shape(type="line",x0=0,y0=y_pos, x1=max_reads + max_reads/20, y1=y_pos,line=dict(dash="dot",color="black",width=2)))
        y_pos=y_pos+len(Annotation_file[Annotation_file["Type"]==type]["ID"])

fig.update_layout(shapes=shapes, barmode='stack',title=go.layout.Title(text="<b>RPKM distribution</b>",xref="paper", x=0.5, y=0.95, font={"size":20}))
fig.update_layout(xaxis=go.layout.XAxis(title=go.layout.xaxis.Title(text="RPKM",
                                                  font=dict(family="Courier New, monospace",size=18,color="#7f7f7f"))))
print("third plot done")
Third_one=plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')

###############################################################################################################################


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

for a in Annotation_file["working"]:
    df=df.append(GeneCount("../05_featureCount/sub_count/"+a +"/FeatureCountTable", a), ignore_index=True)
data=[]
for a in set(Annotation_file["Type"]):
    xx=list(df[df["ID"].str.startswith(a)].mean())
    data.append(go.Scatter(x=["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"], y=xx, name=a))

fig = go.Figure(data=data)
fig.update_layout(title=go.layout.Title(text="<b>Sequence Depth Evaluation</b>",xref="paper", x=0.5, y=0.95, font={"size":20}))
fig.update_layout(xaxis=go.layout.XAxis(title=go.layout.xaxis.Title(text="Fraction of sequenced Reads",
                                                  font=dict(family="Courier New, monospace",size=18,color="#7f7f7f"))))
fig.update_layout(yaxis=go.layout.YAxis(title=go.layout.yaxis.Title(text="Number of Genes > 10 Reads",
                                                  font=dict(family="Courier New, monospace",size=18,color="#7f7f7f"))))

Fourth_one=plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')



##########################################################################################################################
#f = open(Analyses_Folder+'06_Overview_Mapping/Mapping_overview.html','w')
f=open(sys.argv[2], "w")
message = f"""<html>
<h1><u>RNA Pipeline Module 1 Overview</u></h1>
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
<head></head>
<body>{first_one}<br/><hr><br/>{Second_one}<br/><hr><br/>{Third_one}<br/><hr><br/>{Fourth_one}</body>
</html>"""

f.write(message)
f.close()
