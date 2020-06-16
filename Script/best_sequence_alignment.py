def addToSeq(seq,w):
    if w == "begin":
        newSeq= ["-"]+ seq
    elif w == "end":
        newSeq= seq + ["-"]
    return newSeq

def delFromSeq(seq,w):
    if w == "begin":
        del seq[0]
    elif w == "end":
        del seq[-1]
    return seq

def calcScore(matrix,seq1,seq2):
    score=0
    for i in range(len(seq1)):
        pair=seq1[i]+seq2[i]
        if "-" in pair:
            score-=2
        else:
            value=matrix[pair]
            score+=value
    return score

def prepSeq(seq,length,w):
    for i in range (length):
        if len(seq)<length and w == "begin":
            seq=addToSeq(seq,w)
        elif len(seq)<length and w== "end":
            seq=addToSeq(seq,w)
    return list(seq)

def slidingSeq(seq1,seq2,w):
    sequences=[]
    if ("-" not in seq1 and seq2[0]=="-") or ("-" not in seq1 and "-" not in seq2):
        seq2=addToSeq(seq2,"end")
        seq2=delFromSeq(seq2,"begin")
        sequences.append(seq1)
        sequences.append(seq2)
    elif w=="end":
            seq1=delFromSeq(seq1,"end")
            seq2=delFromSeq(seq2,"begin")
            sequences.append(seq1)
            sequences.append(seq2)
    elif w=="begin":
        seq1=addToSeq(seq1,"begin")
        seq2=addToSeq(seq2,"end")
        sequences.append(seq1)
        sequences.append(seq2)
    print(sequences)
    return sequences

def calcAlignment(seq1,seq2,matrix):
    length=len(seq1)+len(seq2)
    if seq1==seq2:
        score=calcScore(matrix,seq1,seq2)
        scoredSequences=[seq1,seq2,score]
        return scoredSequences
    elif len(seq1)<len(seq2):
        seq1=prepSeq(seq1,length,"begin")
        seq2=prepSeq(seq2,length,"end")
    else:
        seq1=prepSeq(seq1,length,"end")
        seq2=prepSeq(seq2,length,"begin")
    score=calcScore(matrix,seq1,seq2)
    scoredSequences=[seq1,seq2,score]
    print(scoredSequences)
    for i in range(len(seq1)):
#            if ("-" not in seq1 and seq2[-1]=="-") or (seq2[-1]=="-"):
            if seq2[-1]=="-":
                w="begin"
                sequences=slidingSeq(seq1,seq2,w)
                seq1=sequences[0]
                seq2=sequences[1]
                newScore=calcScore(matrix,seq1,seq2)
                if newScore>score:
                    scoredSequences[0]=[seq1,seq2,newScore]
                    score=newScore
            else:
                w="end"
                sequences=slidingSeq(seq1,seq2,w)
                seq1=sequences[0]
                seq2=sequences[1]
                newScore=calcScore(matrix,seq1,seq2)
            if newScore>score:
                scoredSequences=[seq1,seq2,newScore]
                score=newScore

    return scoredSequences

seq1=list("ATG")
seq2=list("AC")
matrix={"AA":2,"AT":-1,"AC":-1,"AG":0,
        "TA":-1,"TT":2,"TC":0,"TG":-1,
        "CA":-1,"CT":0,"CC":2,"CG":-1,
        "GA":0,"GT":-1,"GC":-1,"GG":2}

x=calcAlignment(seq1,seq2,matrix)
print(x[0])
print(x[1])
print("score is ", x[2])
