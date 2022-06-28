
file=open("RMVar_Human_basic_info_m6A.txt","r")
data=file.readlines()

def Findmotif(genename):
    fout=open(genename+"_MOTIF.FASTA","w")
    n=1
    for i in data:
        if genename in i:
            j=i.strip().split()
            moti_sequence=j[11]
            seq=moti_sequence[28:37]
            fout.write(">seq"+str(n)+"\n"+str(seq)+"\n")
            n=n+1
        else:
            continue
    fout.close()

list=["LINC01235","P3H2-AS1","LINC01198"]
        
for genename in list:
    Findmotif(genename)
