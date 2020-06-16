#1) build a 0 matrix che ha lunghezza (len(seq1)+1)*(len(seq1)+1)
#2) initialize F(0,0)
#3) fill first column with the gaps penalties (-2,-4,-6, ...)
#4) fill first raw in the same way
#5) for the rest of the position using the maximum F(i,j) -> match, gap1, gap2
seq1="ATCAGGGGATCG"
seq2="ATCAATCG"
matrix={"AA":2,"AT":-1,"AC":-1,"AG":0,
        "TA":-1,"TT":2,"TC":0,"TG":-1,
        "CA":-1,"CT":0,"CC":2,"CG":-1,
        "GA":0,"GT":-1,"GC":-1,"GG":2}
C=len(seq1)+1 #colonne
R=len(seq2)+1 #righe

# Associa ad ogni match tra due nucleo il suo valore della matrice
F=[[0]*C for x in range(R)] # THE VALUE MATRIX
P=[["0"]*C for x in range(R)] # THE COORDINATE MATRIX ("u","l","d")
for i in range(1,len(P[0])):
    P[0][i]="l"
for j in range(1,len(P)):
    P[j][0]="u"
# print(P)
# print(F)

pen=-2 #penalty
for c in range(1, C):
    #print("Col",c)
    for r in range(1, R):
        #print("row",r)
        al=seq1[c-1]+seq2[r-1]
        #print(al)
        diag=F[r-1][c-1]+matrix[al]
        left=F[r][c-1]+pen
        up=F[r-1][c]+pen
        #max1=max(diag,left,up)
        zipp=list(zip((diag,left,up),("d","l","u"))) # The zip here assign to each numerical value of the diagonal, left and up of the F matrix,
                                                     # the letter "d","l" and "u" fot the P matrix... Print it
        # print(zipp)
        # print(r, " ", c)
        F[r][c],P[r][c]=max(zipp) # that allows to find the maximum of the three possibilities of F for each the position... Remember that each of
                                  #of the three value ha its own letter "d","u","l"
        #print(F)
# print(F)
# print(P)

R=R-1
C=C-1
print("the bottom right score is: " ,F[R][C]) #this is the bottom right value, That should be the max score (not always)
#print(P[R][C]) #this is the bottom right coordinate
temp=""
targ=""
n=0
while P[R][C]!="0": # or F[C][R]!=0   -> IDK if it's needed
    n=n+1
    if P[R][C]=="l":
        temp=temp+seq1[C-1]
        targ=targ+"-"
        C=C-1
        # print(temp[::-1])
        # print(targ[::-1])
    elif P[R][C]=="d":
        temp=temp+seq1[C-1]
        targ=targ+seq2[R-1]
        C=C-1
        R=R-1
        # print(temp[::-1])
        # print(targ[::-1])
    elif P[R][C]=="u":
        temp=temp+"-"
        targ=targ+seq2[R-1]
        R=R-1
        # print(temp[::-1])
        # print(targ[::-1])
#print(temp)
print(temp[::-1])
print(targ[::-1])

        #print(al)
        #F[r][c]=matrix[al]
#print(F)
