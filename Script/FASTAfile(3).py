fasta=open("/media/alessandro/DATA/User/BIOINFORMATICA.BOLOGNA/Programming_for_Bioinformatics/Module2/Exercise/sequences.txt","r")
list=[]
l1=[]
n=0
#print(len(l1))
for line in fasta:
    if ">" not in line:
        list.append(line[:-1])
        n+=1
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
    lenght=len(fasta)
    fastlist.append(fasta)
    fastlist.append(lenght)
    fastafin.append(fastlist)
print(fastafin)

min=fastafin[0][1]
max=fastafin[0][1]
sum=0
for el in fastafin:
    sum+=el[1]
    if el[1]<min:
        min=el[1]
        seqmin=el[0]
    elif el[1]>max:
        max=el[1]
        seqmax=el[0]
mean=(float(sum)/len(fastafin))
print("The average lenght of the FASTA sequences is: ", mean)
print("The longest sequence with" ,max, "residues is:",seqmax)
print("The shortest sequence with" ,min, "residues is:",seqmin)
