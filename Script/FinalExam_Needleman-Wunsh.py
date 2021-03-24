import numpy as np #chose to import numpy and use it for the matrices
import input_data1 #I'll import "input_data.py" file where seq1, seq2 and substitution matrix are stored
sub_matrix = input_data1.BLOSUM52 #The substitution matrix that we will use is the BLOSUM52
# print(sub_matrix)
seq1= input_data1.seq1 #I import the seq1
seq2= input_data1.seq2 #I import the seq2
pen=-1 #The penalty value

def mat_creation(sub_matrix,seq1,seq2,pen):
    #a=np.zeros((2,4))
    F=np.zeros((len(seq2)+1,len(seq1)+1)) #Initialization of the scoring matrix, where seq2 define the number of rows in the matrix and seq1 the number of colums
    #F=[[0]*(len(seq1)+1) for x in range(len(seq2)+1)]  # list comprehension
    P=np.zeros((len(seq2)+1,len(seq1)+1),dtype=str) #Initialization of the traceback matrix, where seq2 define the number of rows in the matrix and seq1 the number of columns, and I specified that the values have to be string
    #P=[[""]*(len(seq1)+1) for x in range(len(seq2)+1)] # list comprehension
    for c in range(1,len(P[0])): #filling the first row with "l" that stand for left
        P[0][c]="l"
    for r in range(1,len(P)): #filling the first column with "u" that stand for up
        P[r][0]="u"
    #MATRICES CONSTRUCTION
    #The following loops are for filling up the matrices, from the top left to the buttom row by row (????????????????????????)
    for r in range(1,len(seq2)+1):
        for c in range(1,len(seq1)+1):
            al=seq1[c-1]+seq2[r-1] #The matching pairs for each couple of the two sequences
            diag=F[r-1][c-1]+sub_matrix[al] #The score coming from the diagonal
            left=F[r][c-1]+pen #The score coming from the left, so in the same row but in the previous column
            up=F[r-1][c]+pen #The score coming from the above, so in the same column but in the upper (previous) row
            zipp_list=list(zip((diag,left,up),("d","l","u"))) #This will associate each resolut from diag, letf, up to a string variable "d","l" and "u" in a list that has for each position all the three possible scores and relative derivations
            # print(zipp_list)
            F[r][c],P[r][c]=max(zipp_list) #In this way I assign the maximum value and the relative direction of the three direction saved in the zipp_list, respectively to the F and P matrices
    return(F,P)

def local_align(seq1,seq2,F,P):
    # print(seq1,seq2)
    print(F)
    print(P)
    max=F[0][0] #I initialize the maximum value with the first value of the matrix that will be 0
    list_max=[]
    list_coord=[]
    for r in range(1,len(seq2)+1):
        for c in range(1,len(seq1)+1): #I loop again for every single value of the F matrix to see where is the maximum that will be my starting point
            if max<=F[r][c]:
                max,R,C=F[r][c],r,c
                list_max.append(max)
                list_coord.append([R,C])
                #I assign to max the max value in the matrix, and R and C are the coordinates of the max, but looped on the sequence lenghts
    print(list_max)
    print(list_coord)
    # print(max,R,C) # 3 and 5 are the coordinates in the sequence not in the matrix
    # print(F[R][C])
    temp,targ="","" #this are the two sequence resulting from the alignment
    while F[R][C]!=0: #I keep the iteration till I reach the 0 value in the matrix with the scores
        #print(R,C)
        if P[R][C]=="l": #If in the traceback matrix there's an "l" that means that the score derives from the left cell
            targ=targ+"-"  #That also means that in the targ sequence of the alignent we will have a gap
            temp=temp+seq1[C-1]
            C=C-1 #I will move the "pointer" to the left cell, that is where the score has derived
        elif P[R][C]=="d":  #If is not "l" we check if it comes from the diagonal "d"
            temp=temp+seq1[C-1]
            targ=targ+seq2[R-1]
            C=C-1
            R=R-1 #In this case I have to move one step to the left and un step to the top --> in diagonal
        elif P[R][C]=="u": #I check if it derives from the upper "cell"
            temp=temp+"-"
            targ=targ+seq2[R-1]
            R=R-1 #
    return(max,temp,targ)

F,P = mat_creation(sub_matrix,seq1,seq2,pen)
#print(F)
#print(P)
max,temp,targ=local_align(seq1,seq2,F,P)
print("The best local alignment is: ")
print(temp[::-1])
print(targ[::-1])
print("with a score of :" , max)
