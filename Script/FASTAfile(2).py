import sys

#fasta=open("/media/alessandro/DATA/User/BIOINFORMATICA.BOLOGNA/Programming_for_Bioinformatics/Module2/Exercise/sequences.txt","r")

#print(len(l1))
def read_fasta(fasta):
    list=[]
    l1=[]
    for line in fasta:
        if ">" not in line:
            list.append(line[:-1])
        elif ">" in line:
            list=[]
        if list not in l1:
            l1.append(list)
    #print(l1)
    fastafin=[]
    for li in l1:
        fastlist=[]
        #print(li)
        fasta="".join(li)
        fastlist.append(fasta)
        fastafin.append(fastlist)
    return(fastafin)



if __name__=="__main__":
    fasta1=sys.argv[1]
    fasta=open(fasta1)
    print(read_fasta(fasta))
