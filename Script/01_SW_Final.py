'''This is a program that perform a pairwise local alignment of protein sequences, based on the Smith-Waterman algorithm.
   At first I will import the two sequences as well as the substitution matrix from an external python file. Once imported these:
1) I define the "mat_creation" function with seq1, seq2, substitution matrix and penalty as parameters. This function will firstly
   create the F and P matrices using the numpy package based on the lenght of the sequences. Then the function will loop for each position
   of the F matrix from the top left to the bottom right to fill it up with the best score out of the four possible score of the substring A,B: (r -> row, c -> column)
      1°. F(r-1,c-1) + score substring (Ar, Bc)
      2°. F(r-1,c) + penalty
      3°. F(r,c-1) + penalty
      4°. 0 --> in the Smith Waterman negative scores are not accepted: better to start a new alignment than to have a negative score alignent
   In the same time the function will fill up the P matrix with the position of derivation (up -> "u", left -> "l", diagonal -> "d" or "0") of the maximum score of the
   the substring, this matrix will be used in the next function for the traceback process.
   Both the matrices will be returned to the main program.
2) I define a second function "local_align" with seq1, seq2 and the previously built F and P matrices as parameters. The function will
   find in the F matrix the higher score(s) corresponding to the best alignment(s) and relative coordinate(s) in the matrices. These will be
   our starting point(s) for the while loop, on the F matrix, which will iterate the construction of the aligned sequences ("temp" and "targ" variables) using
   the P matrix for the traceback process, having these 3 cases:
     1° case: there is a "d" in the P matrix, this means that the score in the F(r,c) matrix has derived from the diagonal F(r-1,c-1) and therefore there
              will be a match of the substring
     2° case: there is a "u" in the P matrix, this means that the score has derived from the upper position F(r-1,c) and therefore there is a gap in the seq1 (template sequence)
     3° case: there is a "l" in the P matrix, this means that the score has derived from the left position F(r,c-1) and therefore there is a gap in the seq2 (target sequence)
   With this process I move the r and c indexes following the "instructions" for the directions stored in the P matrix and so the path of the best scoring alignment saved in the
   in the F matrix untill a 0 value in F will be reached, this is the "stop" condition for the Smith-waterman's local alignment process. This while
   loop will perform n times, thanks to a for loop, for each best score found in the F matrix.
   Once build the template and target sequence(s) of the local alignment(s) the function will return them to the main program and the best alignments
   and their score will be printed'''

import numpy as np #chose to import numpy and use it for the generation of the matrices
import input_data #I'll import "input_data.py" file where seq1, seq2 and the substitution matrix are stored

#Import the substitution matrix, the first sequence and the second sequence
sub_matrix = input_data.BLOSUM52 #substitution matrix BLOSUM52
seq1= input_data.seq1 #first sequence, template
seq2= input_data.seq2 #second sequence, target

pen=-1 #define the value of the gap penalty

def mat_creation(sub_matrix,seq1,seq2,pen):
    #MATRICES INITIALIZATION
    F=np.zeros((len(seq2)+1,len(seq1)+1)) #Initialization of the scoring matrix, where seq2 define the number of rows and seq1 the number of columns in the matrix
    P=np.zeros((len(seq2)+1,len(seq1)+1),dtype=str) #Initialization of the traceback matrix (P), where seq2 define the number of rows and seq1 the number of columns, the values of wich it's composef are empty strings
    #It's also possible to genenrate them with the list comprehension method as shown below:
    #F=[[0]*(len(seq1)+1) for x in range(len(seq2)+1)]  # list comprehension
    #P=[[""]*(len(seq1)+1) for x in range(len(seq2)+1)] # list comprehension
    for c in range(1,len(P[0])): #filling up the first row and then the first column of the traceback matrix P, with "0"s
        P[0][c]="0"
    for r in range(1,len(P)):
        P[r][0]="0"

    #MATRICES CONSTRUCTION
    #The following loops are for filling up the matrices, from the top left to the buttom right, row by row with the relative scores
    for r in range(1,len(seq2)+1): #For every row, I iterate from 1 becasue both the first line and first columns should be left just with zeros
        for c in range(1,len(seq1)+1): #For each column/element in the row
            aa_pairs=seq1[c-1]+seq2[r-1] #The matching amino acid pairs of the current posiion in the two sequences (targ and temp). I'll use it for finding the score of the pair in the substitution matrix
            diag=F[r-1][c-1]+sub_matrix[aa_pairs] #The score coming from the diagonal plus the score present in the sub_matrix corresponding to the current aa_pairs
            left=F[r][c-1]+pen #The score coming from the left, so in the same row but in the previous column, sum to the penalty value
            up=F[r-1][c]+pen #The score coming from the above, so in the same column but in the upper row, sum to the penalty value
            zipp_list=list(zip((diag,left,up,0),("d","l","u","0"))) #The zip function will associate the results coming from diag, letf, up and the 0 value respectively to a string variable "d","l","u" and "0" in a list that will have for each position all the four possible scores and relative derivations
            #print("the scores and relative derivations of the cell row:" , r ," column:" ,c, "of the matrix is: ", zipp_list) #To see how is the zipp_list for each position
            F[r][c],P[r][c]=max(zipp_list) #Assign the maximum value, between the four present in the zipp_list and the relative direction of derivation respectively to the F and P matrices, the 0 value will be greater and will take the place of the substring alignment whether this is negative
    return(F,P)

def local_align(seq1,seq2,F,P):

    #DETECTION OF THE BEST ALIGNMENT SCORE(S) (max value(s)) AND ITS COORDINATE(S) for the next alignment procedure
    current_max=F[0][0] #Initialize the maximum value with the first value of the matrix that will be 0, and I will use it for the comparison and the research of the highest alignment score
    list_max=[]
    starting_coord=[]
    for r in range(1,len(seq2)+1):
        for c in range(1,len(seq1)+1): #Loop for every single value of the F matrix for finding the maximum value of the matrix
            if F[r][c]>current_max: #If the current score of the cell in the F matrix is greater than the current_max
                current_max,R,C=F[r][c],r,c #reassign the current_max, the Row and Column coordinate
                list_max[0]=current_max #assign the current_max to the first position of the list_max
                starting_coord[0]=[R,C] #assign the Row and Column coordinate of the current max in the starting_coord list
            elif F[r][c]==current_max: #Check if the current max was already present in tge list_max (if it is equal to the current_max) IDK IF IT'S RIGHT
                current_max,R,C=F[r][c],r,c
                list_max.append(current_max)
                starting_coord.append([R,C])
    print(list_max)
    print(starting_coord)
    max_value=list_max[0]

    #ALIGNMENT PROCESS performed with tracebacking of the P matrix
    final_align=[] #Initialization of the list that will contain all the final sequence alignent(s), template(s) and target(s) that will be returned to the main program
    for li in range(len(starting_coord)): #Iterate the alignment process for each of the staring point coordinate(s) previously detected
        alignment=[]
        R=starting_coord[li][0] #The first coordinate of the sublist will correspond to the Row starting point
        C=starting_coord[li][1] #The second coordinate of the sublist will correspond to the Column starting point
        temp,targ="","" #Initialization of the two sequences (template and target corresponding to seq1 and seq2) that will be filled up with the alignent
        while F[R][C]!=0: #Keep the iteration till in the F matrix a 0 value will be reached, that will be the ending point of the local alignment process
            if P[R][C]=="l": #If in the traceback matrix P there's an "l" that means that the current processing score had derived from the left cell
                targ=targ+"-"  #That also means that in the targ sequence of the alignent there will be a gap
                temp=temp+seq1[C-1] #while in the temp string there will be the corresponding amino acid
                C=C-1 #I will move the "pointer" to the left cell, that is where the score has derived, and this will be the starting point of the next "while" iteration
            elif P[R][C]=="d":  #If is not "l" we check if the score had derived from the diagonal "d"
                temp=temp+seq1[C-1] #Add to both the target and the template the matching pair of amino acid
                targ=targ+seq2[R-1]
                C=C-1
                R=R-1 #In this case I have to move one step to the left and un step to the top, in diagonal
            elif P[R][C]=="u": #I check if it derives from the upper "cell"
                temp=temp+"-" #In this case the template sequence will have a gap
                targ=targ+seq2[R-1] #While the taget the amino acid corresponding to that position
                R=R-1 #I set the coordinate of the next iteration one row above
        alignment.append(temp)
        alignment.append(targ) #I save both the aligned sequences (template and target) into a "temporary" list
        final_align.append(alignment) #save the above list(s) with the sequences to the final_align list, that will be returned to the main program
    return(max_value,final_align)

#FUNCTION RECALL
F,P = mat_creation(sub_matrix,seq1,seq2,pen)
#print(F)
#print(P)
max_score,final_align = local_align(seq1,seq2,F,P)
#print(final_align, max_score)

#OUTPUT
for el in range(len(final_align)): #Iterate for each alignment saved in the final_align list and print it with relative score(s)
    print("The best" , el+1, "local alignment with a score of", max_score, "is: ")
    print("Temp:",final_align[el][0][::-1])
    print("Targ:",final_align[el][1][::-1])
