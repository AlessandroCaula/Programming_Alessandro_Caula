matrix=open("/media/alessandro/DATA/User/BIOINFORMATICA.BOLOGNA/Programming_for_Bioinformatics/Module2/Exercise/PAM250.txt","r")
matrix1=open("/media/alessandro/DATA/User/BIOINFORMATICA.BOLOGNA/Programming_for_Bioinformatics/Module2/Exercise/PAM250(1).txt","r")
rows="ARNDCQEGHILKMFPSTWYV"
cols="ARNDCQEGHILKMFPSTWYV"
def tr_matrx(matrix,rows,cols):
    dict={}
    scores=[]
    for row in matrix:
        score=[]
        row.rstrip()
        score=row.split()
        scores.append(score)
    #print(scores)
    row=0
    for r in scores:
        col=0
        for c in r:
            k=rows[row]+cols[col]
            dict[k]=float(scores[row][col])
            col+=1
        row+=1
    return(dict)
print(tr_matrx(matrix,rows,cols))

def dic_matrix(file):#funzione che crea il dizionario a partire da matrice di sostituzione
    aminoacid = file.readline().split()
    ami = aminoacid[3] #indicizzato a tre perche' se guardi il file in quella pos della lista della prima riga ci sono gli amminoacidi
    #print(ami)
    couple = {}
    for a in ami:
        col = file.readline().split()#scorro le righe a partire dalla prima con un solo elemento print col
        #print(col)
        for i in range(len(col)):
            couple[a + ami[i]] = float(col[i])#posso anche mettere int(col[i][:-1]) che mi consente di rimuovere il punto
            #che altrimenti da errore se non lo tolgo mettendo int
    file.close()
    return couple
print(dic_matrix(matrix1))
