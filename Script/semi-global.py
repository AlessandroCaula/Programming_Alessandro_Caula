import numpy as np #chose to import numpy and use it for the generation of the 
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
            zipp_list=list(zip((diag,left,up),("d","l","u"))) #The zip function will associate the results coming from diag, letf, up and the 0 value respectively to a string variable "d","l","u" and "0" in a list that will have for each position all the four possible scores and relative derivations
            #print("the scores and relative derivations of the cell row:" , r ," column:" ,c, "of the matrix is: ", zipp_list) #To see how is the zipp_list for each position
            F[r][c],P[r][c]=max(zipp_list) #Assign the maximum value, between the four present in the zipp_list and the relative direction of derivation respectively to the F and P matrices, the 0 value will be greater and will take the place of the substring alignment whether this is negative
    return(F,P)

def local_align(seq1,seq2,F,P):

    #DETECTION OF THE BEST ALIGNMENT SCORE(S) (max value(s)) AND ITS COORDINATE(S) for the next alignment procedure
    current_max=F[0][0] #Initialize the maximum value with the first value of the matrix that will be 0, and I will use it for the comparison and the research of the highest alignment score
    list_max=[]
    list_coord=[]
    for r in range(1,len(seq2)+1):
        for c in range(1,len(seq1)+1): #Loop for every single value of the F matrix for finding the maximum value of the matrix
            if current_max<=F[r][c]: #If the current score of the cell in the F matrix is equal or greater than the current_max
                current_max,R,C=F[r][c],r,c #reassign the current_max, the Row and Column coordinate
                list_max.append(current_max) #save and store the current_max in the list_max list
                list_coord.append([R,C]) #save and store the Row and Column coordinate of the current max in the list_coord list
    #In this way I will end up with all the scores (and their coordinates) in the matrix that during the control resulted equal or greater than the previous current_max. Next step will be to select only the actual maximum(s)
    max_value=max(list_max) #Finding the best score alignment present in the list_max list
    starting_coord=[] #this will be the final list with only the coordinate(s) of the actual maximum score(s), that will be the starting point(s) for the alignment process(es)
    for i in range(len(list_max)):
        if max_value==list_max[i]: #Using the max_value that I've previously exctracted from the same list_max I check if there is another value equal to this apart from the known one, if yes that will be my second starting point
            starting_coord.append(list_coord[i]) #Appending the coordinate of the max vale(s) to this starting_coord list

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
