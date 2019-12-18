seq1="AAATTA"
seq2="TTTTTT"
def score(seq1,seq2):
    sum=0
    if len(seq1)==len(seq2):
        for i in range(len(seq1)):
            if seq1[i]==seq2[i]:
                sum=sum+1
            else:
                sum=sum-1
        print(sum)
    else:
        print("the two sequences have different lenght")
score(seq1,seq2)
