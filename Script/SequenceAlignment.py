file=open("/media/alessandro/DATA/User/BIOINFORMATICA.BOLOGNA/Programming for Bioinformatics/Module2/Exercise/Blosum.table","r")
def table(file):
    lst=[]
    a=file.readline()
    #listc=a.split()
    for line in file:
        c=line.split()
        lst.append(c[1:])
    return(lst)
    file.close()

file1=open("/media/alessandro/DATA/User/BIOINFORMATICA.BOLOGNA/Programming for Bioinformatics/Module2/Exercise/Blosum.table","r")
def dictionary(lst1,file1):
    dict={}
    amino=file1.readline()
    #print(amino)
    listc=amino.split()
    #print(listc)
    for i in range(len(listc)):
        for j in range(len(listc)):
            c=listc[i]+listc[j]
            dict[c]=lst1[i][j]
    return(dict)

seq=open("/home/alessandro/Desktop/sequence.txt","r")
def sequences(seq):
    seqtot=[]
    for i in seq:
        #seqli=[]
        if i[0]!= ">":
            line=i.rstrip()
            #seqli.append(line)
            seqtot.append(line)
    return(seqtot)
    seq.close()

def score(dictio,seque):
    score=0
    for i in range(len(seque[0])):
        str=seque[0][i]+seque[1][i]
        #print(str)
        if "-" not in str and str!="--":
            num=int(dictio[str])
            score=num+score
        if "-" in str and str!="--":
            score=score-2
    print(seque[0])
    print(seque[1])
    print(score)


seque=sequences(seq)
print(seque)
x=table(file)
# print(x)
dictio=dictionary(x,file1)
print(dictio)
score(dictio,seque)
