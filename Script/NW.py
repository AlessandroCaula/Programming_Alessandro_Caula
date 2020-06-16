file=open("/media/alessandro/DATA/User/BIOINFORMATICA.BOLOGNA/Programming_for_Bioinformatics/Module2/Exercise/Blosum_table.txt","r")
seq1="MAGPATQS" #PMKLMALQLLLWHSALWTVQEATPLGPASSLPQSFLLKCLEQVRK
seq2="MD" #VIILLQYFLFASLVQSAPVGLPSNLPPAFRETAVEAKTLGDKILKDIS

def dictionary(file):
    dict={}
    listc=[]
    lst=[]
    a=file.readline()
    listc=a.split()
    for line in file:
        c=line.split()
        lst.append(c[1:])
    for i in range(len(listc)):
        for j in range(len(listc)):
            c=listc[i]+listc[j]
            dict[c]=int(lst[i][j])
    return(dict)

C=len(seq1)+1 # Columns
R=len(seq2)+1 # Rows

def mat_score(C,R):
    F=[[0]*C for x in range(R)] # THE VALUE MATRIX
    return(F)

def mat_dir(C,R):
    P=[["0"]*C for x in range(R)] # THE COORDINATE MATRIX ("u","l","d")
    for i in range(1,len(P[0])):
        P[0][i]="l"
    for j in range(1,len(P)):
        P[j][0]="u"
    return(P)

dictionarys=dictionary(file)
#print(dictionarys)
mat_scores=mat_score(C,R)
#print(mat_scores)
mat_dirs=mat_dir(C,R)
#print(mat_dirs)

# MATRICES CONSTRUCTION
pen=-2 #penalty
for c in range(1, C):
    for r in range(1, R):
        al=seq1[c-1]+seq2[r-1]
        diag=mat_scores[r-1][c-1]+dictionarys[al]
        #print(diag)
        left=mat_scores[r][c-1]+pen
        up=mat_scores[r-1][c]+pen
        zipp=list(zip((diag,left,up),("d","l","u"))) # The zip here assign to each numerical value of the diagonal, left and up of the F matrix,                                              # the letter "d","l" and "u" fot the P matrix... Print it
        mat_scores[r][c],mat_dirs[r][c]=max(zipp) # that allows to find the maximum of the three possibilities of F for each the position... Remember that each of
                                  #of the three value ha its own letter "d","u","l"
print(mat_scores)
print(mat_dirs)

R=R-1
C=C-1
print("the bottom right score is: " ,mat_scores[R][C]) #this is the bottom right value, That should be the max score (not always)
temp=""
targ=""
n=0
while mat_dirs[R][C]!="0": # or F[C][R]!=0   -> IDK if it's needed
    if mat_dirs[R][C]=="l":
        temp=temp+seq1[C-1]
        targ=targ+"-"
        C=C-1
    elif mat_dirs[R][C]=="d":
        temp=temp+seq1[C-1]
        targ=targ+seq2[R-1]
        C=C-1
        R=R-1
    elif mat_dirs[R][C]=="u":
        temp=temp+"-"
        targ=targ+seq2[R-1]
        R=R-1
print(temp[::-1])
print(targ[::-1])
