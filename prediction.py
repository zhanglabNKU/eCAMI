#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 20:08:33 2019

@author: joy1314bubian
"""

import numpy as np

import os
import datetime
import re
import psutil,argparse,sys
from multiprocessing import Pool 



def read_file(database_dir,input_fasta_file):
    f = open(database_dir +input_fasta_file, 'r')
    text=f.read()
    f.close()
    text=text.split('\n')
    while '' in text:
        text.remove('')
    proteinname=[]
    new_text={}

    for i in range(len(text)):
        if '>' in text[i]:
            proteinname.append(text[i])
        else:
            if '>' not in text[i] and '>' in text[i-1]:
                new_text[proteinname[len(proteinname)-1]]=text[i]
            else:
                new_text[proteinname[len(proteinname)-1]]+=text[i]
    protein_string=[]
#    print(len(proteinname))
    for each_name in proteinname:
        protein_string.append(new_text[each_name])
    return [proteinname,protein_string]


def find_all(sub,s):
	index_list = []
	index = s.find(sub)
	while index != -1:
		index_list.append(index+1)
		index = s.find(sub,index+1)
	
	if len(index_list) > 0:
		return index_list
	else:
		return -1



def get_cluster_number(file_name,n_mer_dir_path,output_dir,n_mer,piece_number,protein_name,protein_string,important_n_mer_number,beta):

    kmer_dir_list=os.listdir(n_mer_dir_path)

    for i in range(len(kmer_dir_list)):
        kmer_dir_list[i]=kmer_dir_list[i]+'/kmer_for_each_cluster'
    dir_name=n_mer_dir_path

    kmer_label=[]
    kmer_message=[]
    all_cluster_kmer=[]
    for each_kmer_dir in kmer_dir_list:
        kmer_file_list=os.listdir(dir_name+'/'+each_kmer_dir)
        for each_file in kmer_file_list:
            f=open(dir_name+'/'+each_kmer_dir+'/'+each_file,'r')
            text=f.read()
            text=text.split('\n')
            f.close()
            kmer_label.append(each_kmer_dir.split('/')[0]+':'+text[0].split(':')[1])

            kmer_message.append(text[1].split('\t')[1])
            each_cluster_kmer=text[2]
            each_cluster_kmer=each_cluster_kmer.split('\t')
            while '' in each_cluster_kmer:
                each_cluster_kmer.remove('')
            cluster_kmer={}
            for each_kmer in each_cluster_kmer:
                cluster_kmer[each_kmer.split(':')[0]]=float(each_kmer.split(':')[1])
            all_cluster_kmer.append(cluster_kmer)
    f=open(output_dir+file_name,'w')
    for i in piece_number:
#        print(i)
        kmer=[]
        for j in range(n_mer):
            text=protein_string[i][j:len(protein_string[i])]
            kmer+=re.findall('.{'+str(n_mer)+'}', text)
        kmer=list(set(kmer))
        each_cluster_score=[]
        each_cluster_number_score=[]
        each_same_kmer=[]
        for j in range(len(all_cluster_kmer)):
            score=0
            number_score=0
            each_cluster_kmer=all_cluster_kmer[j]
            same_kmer=list(set(kmer).intersection(set(each_cluster_kmer.keys())))
            score+=len(same_kmer)
            if score>=important_n_mer_number:
                for k in same_kmer:
                    number_score+=each_cluster_kmer[k]
            each_cluster_score.append(score)
            each_cluster_number_score.append(number_score)
            each_same_kmer.append(same_kmer)
        max_similary_cluster_label=np.argmax(each_cluster_score)
#        py = psutil.Process(os.getpid())
#        print(str(round(py.memory_info()[0]/2.**30,2))+':GB')

        if each_cluster_number_score[max_similary_cluster_label]>=beta and each_cluster_score[max_similary_cluster_label]>=important_n_mer_number:

            f.write(protein_name[i]+'\t'+str(kmer_label[max_similary_cluster_label])+'\t'+str(kmer_message[max_similary_cluster_label])+'\n')
            kmer_index_dict={}
            for each_kmer in each_same_kmer[max_similary_cluster_label]:
                index=find_all(each_kmer,protein_string[i])
                kmer_index_dict[np.min(index)]=[each_kmer,index]
            sorted_kmer_index_dict=sorted(kmer_index_dict.items(),key=lambda x:x[0])
            for each_member in sorted_kmer_index_dict:
                f.write(each_member[1][0])
                f.write('(')
                index=each_member[1][1]
                for m in range(len(index)):
                    if m==len(index)-1:
                        f.write(str(index[m]))
                    else:
                        f.write(str(index[m])+',')
                f.write(')'+',')
            f.write('\n')
    f.close()






def get_validation_results(input_fasta_file,database_dir,output_dir,output_file_name,n_mer,important_n_mer_number,beta,n_mer_dir_path,jobs):
    if database_dir:
        database_dir=database_dir+'/'
    [protein_name,protein_string] = read_file(database_dir,input_fasta_file)

    current_path=os.getcwd()

    all_output_dir_name=output_dir.split('/')
    while '' in all_output_dir_name:
        all_output_dir_name.remove('')
    for i in range(len(all_output_dir_name)):
        if i==0:
            temp_dir_list=os.listdir(current_path)
            if all_output_dir_name[i] not in temp_dir_list:
                os.mkdir(all_output_dir_name[i])
        else:
            temp_dir_name='/'.join(all_output_dir_name[0:i])
            temp_dir_list=os.listdir(temp_dir_name)
            if all_output_dir_name[i] not in temp_dir_list:
                os.mkdir(temp_dir_name+'/'+all_output_dir_name[i])
    temp_output_dir='/'.join(all_output_dir_name)  
    pred_file_name=output_file_name
    if temp_output_dir:
        temp_output_dir=output_dir+'/'

    pool = Pool(processes=jobs)
    piece_number=np.array_split(list(range(len(protein_name))),jobs)

    for i in range(jobs):
        file_name=input_fasta_file+'_thread_'+str(i)+'.txt'
#        get_cluster_number(file_name,n_mer_dir_path,n_mer_dir_name,temp_output_dir,n_mer,piece_number[i],protein_name,protein_string,important_n_mer_number,beta)

        pool.apply_async(get_cluster_number, (file_name,n_mer_dir_path,temp_output_dir,n_mer,piece_number[i],protein_name,protein_string,important_n_mer_number,beta,))
    pool.close()
    pool.join()
    fw=open(temp_output_dir+pred_file_name,'w')
    for i in range(jobs):
        f=open(temp_output_dir+input_fasta_file+'_thread_'+str(i)+'.txt')
        text=f.read()
        f.close()
        fw.write(text)
        os.remove(temp_output_dir+input_fasta_file+'_thread_'+str(i)+'.txt')
    fw.close()
        






def arg_parser(args): 
    '''
    Process command-line user arguments and prepare for further directories
    '''

    ####################################
    #Folders and files:
    parser = argparse.ArgumentParser(description='The input parameters for the identification technology')

    parser.add_argument('-kmer_db',     default="CAZyme",type=str,        help="Change n_mer directories path for prediction")
    parser.add_argument('-output',     default="examples/prediction/output/test_pred_cluster_labels.txt",type=str,        help="file name for prediction results saving")
    parser.add_argument('-input',     default="examples/prediction/input/test.faa",type=str,        help="Define the fasta file name")
    parser.add_argument('-k_mer',           default=8,         type=int,                 help="Peptide length for prediction")
    parser.add_argument('-jobs',            default=8,         type=int,                 help='Number of processor for use for prediction')
    parser.add_argument('-important_k_mer_number',        default=3,      type=int,               help="Minimum number of n_mer for prediction")
    parser.add_argument('-beta',        default=1,      type=float,               help="Minimum sum of percentage of frequency of n_mer for prediction")


    args = parser.parse_args(args)

    return (args)







            
if __name__ == "__main__":
    starttime = datetime.datetime.now()
    args = arg_parser(sys.argv[1:])

    inputs=args.input
    inputs=inputs.split('/')
    database_dir='/'.join(inputs[0:len(inputs)-1])
    input_fasta_file=inputs[-1]
    outputs=args.output
    outputs=outputs.split('/')
    
    output_dir='/'.join(outputs[0:len(outputs)-1])
    output_file_name=outputs[-1]

    n_mer=args.k_mer
    n_mer_dir_path=args.kmer_db

    n_mer_dir_path=n_mer_dir_path.split('/')
    while '' in n_mer_dir_path:
        n_mer_dir_path.remove('')
    n_mer_dir_path='/'.join(n_mer_dir_path)    


    important_n_mer_number=args.important_k_mer_number
    jobs=args.jobs
    beta=args.beta

    print(input_fasta_file)
    get_validation_results(input_fasta_file,database_dir,output_dir,output_file_name,n_mer,important_n_mer_number,beta,n_mer_dir_path,jobs)

    endtime = datetime.datetime.now()
    print('total time:'+str((endtime - starttime).total_seconds())+'s')
