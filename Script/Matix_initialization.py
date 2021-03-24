# def score_m(seq1,seq2,ScorongMAtrix,Penalties)
# intialize F to the (0)null matrix
# ini F(i,0)
# ini F(0,j)
# itaration
# return F
seq1="ATCCT" #numer of columns
seq2="ATCT" #number of rows

tot_mat=[]
rows=1
for row in range(1,len(seq2)+1): # seq 2 will define the number of rows, and the number of substring
    mat=[] # For each row I will create a new substring
    for col in range(1,len(seq1)+1): # seq 1 instead will define the columns, and the number 0s (or whatever) in the sbustring
        mat.append(rows) # Append for each column the 0 (or whatever)
    tot_mat.append(mat) # Will append each new rown, already full of 0s columns. This is iterated from the first for cycle, defined by the lenght of seq2
    rows+=1
print(tot_mat)


M=len(seq1)+1
N=len(seq2)+1
F=[[0]*M for x in range(N)]
print(F)

def zero_matrix(seq1,seq2):
    null_matrix=[]
    for i in range(len(seq2)+1):
        row_matrix=[]
        for j in range(len(seq1)+1):
            row_matrix.append(0)
        null_matrix.append(row_matrix)
    return(null_matrix)

def init_m(matrix):
    for i in range(1,len(matrix[0])):
        matrix[0][i]=matrix[0][i-1]-2
    for j in range(1,len(matrix)):
        for n in range(len(matrix[j])):
            matrix[j][0]=matrix[j-1][0]-2
    return(matrix)

# def score(seq1,seq2):
#     for i in range(len(seq1)):
#         if seq1[i]==seq2[i]:



matrix=zero_matrix(seq1,seq2)
print(matrix)
matrix1=init_m(matrix)
print(matrix1)
