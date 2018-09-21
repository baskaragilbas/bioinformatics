
#Gilang Baskara Putera Anyra - 16/394084/PA/17175

import pandas 
import numpy 

#load nilai blosum
blosum = pandas.read_csv('blosum62.csv')


# Fungsi needleman
def needleman(a,b):
    row = len(b)+1
    col = len(a)+1
    #Matrix val = value, destination = arah, blos - Blosum Score.
    val = numpy.zeros(shape=(row,col), dtype=numpy.int)
    destination = numpy.full(shape=(row,col), fill_value=" ", dtype=numpy.str)
    blos = numpy.zeros(shape=(row,col), dtype=numpy.int)
    #score gap
    gap = 6
    #pengisian matrix 
    val[0][0]=0
    for i in range(1,row):
        val[i][0] = val[i-1][0]-gap
    for j in range(1,col):
        val[0][j] = val[0][j-1]-gap    
    for i in range (1,row):
        for j in range (1,col):
            match = val[i-1][j-1]+blosum[a[j-1]][b[i-1]]
            gap_a = val[i][j-1]-gap
            gap_b = val[i-1][j]-gap
            #memilih value maksimum dari kemungkinan yang ada, dan diberi arah
            val[i][j] = max(match,gap_a,gap_b) 
            if max(match,gap_a,gap_b) ==match:
                destination[i][j] = '↖'
            elif max(match,gap_a,gap_b) == gap_a:
                destination[i][j] ='←'
            elif max (match,gap_a,gap_b)==gap_b:
                destination[i][j] ='↑'
    #inisialisasi matrix arah
    destination[0][0] = "-"
    for i in range(1, row):
        destination[i][0] = "↑"
    for j in range(1, col):
        destination[0][j] = "←"
    #pengisian matrix blos dengan blosumscore
    blos[0][0]=0
    for i in range(1,row):
        blos[i][0] = blos[i-1][0]-gap
    for j in range(1,col):
        blos[0][j] = blos[0][j-1]-gap 
    for j in range(1,col):
        for i in range(1,row):
            blos[i][j] = blosum[a[j-1]][b[i-1]]
    #traceback
    align_a = ''
    align_b = ''
    i = len(b)
    j = len(a)
    score = 0
    while i>0 or j>0:
        if destination[i][j] == '↖':
            align_a = a[j-1]+align_a
            align_b = b[i-1]+align_b
            score=score+blos[i][j]
            i = i-1
            j = j-1         
        elif destination[i][j] == '←':
            align_a = a[j-1]+align_a
            align_b = '_'+align_b
            score=score-gap
            j = j-1
        elif destination[i][j] == '↑':
            align_a = '_'+align_a
            align_b = b[i-1]+align_b
            score=score-gap
            i=i-1
        
        elif j>0:
                while j>0:
                    align_a = a[j-1]+align_a
                    align_b = '_'+align_b
                    score=score-gap
                    j=j-1
        elif i>0:
                while i>0:
                    align_a = '_'+align_a
                    align_b = b[i-1]+align_b
                    score=score-gap
                    i=i-1
        else:
            break
        
    print(align_a)
    print(align_b)
    print("SCORE:",score)
    print("\n")
    print(val)
    print("\n")
    print(destination)
    print("\n")
    print(blos)

#Testing
needleman('MNALQM','NALMSQA')


