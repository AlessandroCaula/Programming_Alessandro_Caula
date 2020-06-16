file=open("/media/alessandro/DATA/User/BIOINFORMATICA.BOLOGNA/Programming_for_Bioinformatics/Module2/Exercise/Blosum_table.txt","r")
def dictionary(file):
    dict={}
    listc=[]
    lst=[]
    a=file.readline()
    listc=a.split()
    for line in file:
        c=line.split()
        lst.append(c[1:])
    #print(lst)
    for i in range(len(listc)):
        for j in range(len(listc)):
            c=listc[i]+listc[j]
            dict[c]=lst[i][j]
    return(dict)

dictionary=dictionary(file)
print(dictionary)
