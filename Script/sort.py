#ORDINARY SORT
A=[7,6,5,4,3,2,1]
for j in range(1,len(A)):
    key=A[j]
    i=j-1
    #j=j+1
    while i>-1 and A[i]>key:
        A[i+1]=A[i]
        #print(A[i+1], A[i])
        i=i-1
        #print(i)
    #print(i)
    A[i+1]=key
print(A)

#MERGE SORT
A=[3,2,4,5,6,1]
# p=0
# r=len(A)
# def split(A,p,r):
#     # if p>r:
#     #     return(A)
#     # else:
#     if p<r:
#         q=[len(A)/2]
#         #print(q)
#         return(split(A,p,q),split(A,q+1,r))
#         #return(split(A,q+1,r))
# split(A,p,r,)
p=0
r=len(A)
def split(A):
    if len(A)==1:
        ris=A[0]
        print(ris)
    else:
        half=len(A)//2
        first=A[:half]
        second=A[half:]
        split(first)
        split(second)
print(split(A))
