# -*- coding: utf-8 -*-
"""
@author: Minjeong Kang
"""
# Description: rocesses variant mapping data from multiple files, consolidates mutation information with associated values, and writes the results to a CSV file.
# 설명: 여러 파일에서 변이 매핑 데이터를 처리하고, 변이 정보와 관련 값을 통합하여 CSV 파일로 저장.

WT2_dic = {}
mutation = []
for i in range(0,17):
    a=i*100
    f=open('WG2_' + str(a) + '_Lib_1 trimmed mapping variants.csv', 'r')
    f.readline()
    temp=f.readlines()
    temp2=[]
    for i in temp:
        j=i.replace('"','')
        temp2.append(j.split(','))    
    for x in temp2:
        WT2_dic[(a,str(x[0]+x[1]+x[2]+x[3]+x[4]))]= [x[9], x[13:]]
        mutation.append(str(x[0]+x[1]+x[2]+x[3]+x[4]))
    f.close()
    
for i in range(1,17):
    a=i*100
    f=open('WG2_' + str(a) + '_Lib_1 trimmed mapping variants.csv', 'r')
    f.readline()
    temp=f.readlines()
    temp2=[]
    for i in temp:
        j=i.replace('"','')
        temp2.append(j.split(','))    
    for x in temp2:
        
        WT2_dic[(a+800,str(x[0]+x[1]+x[2]+x[3]+x[4]))]= [x[9], x[13:]]
        mutation.append(str(x[0]+x[1]+x[2]+x[3]+x[4]))
    f.close()
    
mutation = list(set(mutation))

f= open('WT2 variants.csv', 'w')
for i in mutation:
    f.write(i)
    for j in range (0,17):
        a = j*100
        try:
            f.write(','+ WT2_dic[(a,i)][0])
            gene = str(WT2_dic[(a,i)][1])
        except:
            f.write(','+ '0')
    for k in range (1,17):
        b = k*100+800
        try:            
            f.write(','+ WT2_dic[(b,i)][0])
            gene = str(WT2_dic[(b,i)][1])
        except:
            f.write(','+ '0')
    f.write(','+ gene + '\n')
f.close()

    
            