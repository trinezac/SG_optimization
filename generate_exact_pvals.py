# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 21:09:56 2023

@author: Anders
"""
import math

def S2(n,k):
    return sum((-1)**i * math.comb(k, i) * (k-i)**n for i in range(k+1))//math.factorial(k)

def CCP_pval_gen(n,k,d):
    pval=(math.comb(n,d)*S2(k,d)*math.factorial(d))/(n**k)
    
    return(pval)


S2_number_dir="C:/Users/Anders/Downloads/S2_numbers_ordered.tsv"

out="C:/Users/Anders/Downloads/pvals_n_k_d.tsv"


with open(out,'w') as wh:
    wh.write('n\tk\td\tpval\n')
    for n in range(90,101):
        print(n)
        for k in range(0,3000):
            #the chance of drawing 101 balls different balls out of 100 balls will be zero, so not worth spending time
            #testing here.
            upper_lim=min([k,n])
            newline=""
            for d in range(0,upper_lim+1):
                pval=CCP_pval_gen(n,k,d)
                
                newline=newline+str(n)+'\t'+str(k)+"\t"+str(d)+"\t"+str(round(pval,100))+'\n'
            wh.write(newline)
