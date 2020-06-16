#I open the file, and loop it to read each line, until if 1° element == to ATOM, Then I loop it to make the list of the lists with the distances
#RSMD formula --> sqrt((1/n)*sommatoria)
file=open("/media/alessandro/DATA/User/BIOINFORMATICA.BOLOGNA/Programming for Bioinformatics/Module2/Exercise/model8.pdb","r")
file1=open("/media/alessandro/DATA/User/BIOINFORMATICA.BOLOGNA/Programming for Bioinformatics/Module2/Exercise/model8.1.pdb","r")
from math import *
def read_line(file):
    list=[]
    for line in file:
        if "ATOM" and "CA" in line:
            coord=line[32:55]
            coord1=coord.split()
            list.append(coord1)
    return(list)
list1=read_line(file)
list2=read_line(file1)
print(list1)
print(list2)

def RMSD(list1,list2):
    rmsd=0
    for i in range(len(list1)):
        rmsd=rmsd+(float(list1[i][0])-float(list2[i][0]))**2+(float(list1[i][1])-float(list2[i][1]))**2+(float(list1[i][2])-float(list2[i][2]))**2
    sq=sqrt(1/len(list1) * rmsd)
    print(sq)
print(RMSD(list1,list2))
