
import matplotlib.pyplot as plt
import math
import statistics
# import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
# from matplotlib.ticker import PercentFormatter

# f=open('GSE79578_bulkRNAseq.txt')
# new_f=open('GSE79578_bulkRNAseq_96h.bed','a')
# # for line in f.readlines(3):
# #     l=line.strip().split()
# # #     for i,j in enumerate(l):
# # #         print(i,j)
# for line in f.readlines():
#     L = line.strip().split()
# #     # print(L[0],L[10])
#     new_f.write(L[0]+'\t'+L[10]+'\n')

def adding_mathyalation_to_genes_peek_file():
    genes_peek_file=open('genes_peeks_4d_RA_H3K36me3_Intersect.ChIP.peaks.CpG.bed_.bed')
    new_peek_genes_file=open('genes_peeks_4d_RA_H3K36me3_Intersect.ChIP.peaks.CpG.bed','a')
    methyl_file=open('new_4d_RA_H3K36me3Intersect.ChIP.peaks.CpG.bed')

    # genes_peek_file = open('genes_peeks_4d_RA_H3K36me3_Intersect.WCE.ChIP.peaks.CpG.bed_.bed')
    # new_peek_genes_file = open('genes_peeks_4d_RA_H3K36me3_Intersect.WCE.ChIP.peaks.CpG.bed', 'a')
    # methyl_file = open('new_4d_RA_H3K36me3Intersect.WCE.ChIP.peaks.CpG.bed')

    # genes_peek_file = open('genes_peeks_no_2i_2d_H3K36me3_Bismark_PE_Intersect.WCE.ChIP.peaks.CpG.bed_.bed')
    # new_peek_genes_file = open('genes_peeks_no_2i_2d_H3K36me3_Bismark_PE_Intersect.WCE.ChIP.peaks.CpG.bed', 'a')
    # methyl_file = open('new_no_2i_2d_H3K36me3_Bismark_PEIntersect.WCE.ChIP.peaks.CpG.bed')
    #
    # genes_peek_file = open('genes_peeks_no_2i_2d_H3K36me3_Bismark_PE_Intersect.ChIP.peaks.CpG.bed_.bed')
    # new_peek_genes_file = open('genes_peeks_no_2i_2d_H3K36me3_Bismark_PE_Intersect.ChIP.peaks.CpG.bed', 'a')
    # methyl_file = open('new_no_2i_2d_H3K36me3_Bismark_PEIntersect.ChIP.peaks.CpG.bed')


    methyl_percent={}
    peek_methyl_percent={}
    peeks={}
    for line in methyl_file:
        L = line.strip().split()
        peek=L[0]
        methyl=L[1]
        non_methyl=L[2]
        reads=L[3]
        percent=L[4]
        methyl_percent[peek]=(peek,methyl,non_methyl,reads,percent)
    for line_peek in genes_peek_file.readlines():
        L_peek = line_peek.strip().split()
        peek=L_peek[2]
        peeks[peek]=(L_peek[0],L_peek[1],L_peek[2],L_peek[3],L_peek[4],L_peek[5])
    # print(peeks)
    for peek,i in peeks.items():

        for peek_percent_m,j in methyl_percent.items():
            if peek==peek_percent_m:
                if peek not in peek_methyl_percent:
                    peek_methyl_percent[peek]=(i[0],i[1],i[2],i[3],i[4],i[5],j[1],j[2],j[3],j[4])
    for peek,v in peek_methyl_percent.items():
        new_peek_genes_file.write(v[0]+'\t'+v[1]+'\t'+v[2]+'\t'+v[3]+'\t'+v[4]+'\t'+v[5]+'\t'+v[6]+'\t'+v[7]+'\t'+v[8]+'\t'+v[9]+'\n')
    # print(peek_methyl_percent)
# print(adding_mathyalation_to_genes_peek_file())

def met_vs_exp(): ###################################################
    # methylation_file = open('genes_peeks_4d_RA_H3K36me3_Intersect.ChIP.peaks.CpG.bed')
    # gene_exp_file = open('GSE79578_bulkRNAseq_96h.bed')
    # methyl_vs_exp = open('genes_peeks_4d_RA_H3K36me3_Intersect.ChIP.peaks.CpG.metVSexp.bed', 'a')

    methylation_file = open('genes_peeks_4d_RA_H3K36me3_Intersect.WCE.ChIP.peaks.CpG.bed')
    gene_exp_file = open('GSE79578_bulkRNAseq_96h.bed')
    methyl_vs_exp = open('genes_peeks_4d_RA_H3K36me3_Intersect.WCE.ChIP.peaks.CpG.metVSexp.bed', 'a')


    met_dict = {}
    exp_dict = {}
    for line in methylation_file.readlines():
        L = line.strip().split()
        gene_name=L[1]
        met=float(L[6])
        non_met=float(L[7])
        if gene_name not in met_dict:
            met_dict[gene_name]=(met,non_met)
        else:
            methyl,non_methyl=met_dict[gene_name]
            met_dict[gene_name]=(float(methyl+met),float(non_methyl+non_met))
    for exp_line in gene_exp_file.readlines():
        exp_L= exp_line.strip().split()
        gene_name_exp=exp_L[0]
        gene_exp=exp_L[1]
        exp_dict[gene_name_exp]=gene_exp
    for gene_name in met_dict:
        me,non=met_dict[gene_name]
        if gene_name in exp_dict.keys():
            methyl_vs_exp.write(gene_name+'\t'+exp_dict[gene_name]+'\t'+str(me)+'\t'+str(non)+'\t'+(str(float(me)/(float(me+non))))+'\n')
# print(met_vs_exp())

def chip_nimus_wce():
    chip_methyl_vs_exp_file = open('genes_peeks_4d_RA_H3K36me3_Intersect.ChIP.peaks.CpG.metVSexp.bed')
    wce_methyl_vs_exp_file = open('genes_peeks_4d_RA_H3K36me3_Intersect.WCE.ChIP.peaks.CpG.metVSexp.bed')
    files=[chip_methyl_vs_exp_file,wce_methyl_vs_exp_file]
    chip_nimus_wce_file=open('genes_peeks_4d_RA_H3K36me3_Intersect.WCE.ChIP.peaks.CpG.metVSexp_delta.bed','a')
    chip_nimus_wce_file.write(
        'gene' + '\t' + 'exp_wce' + '\t' + 'met_wce' + '\t' + 'non_met_wce' + '\t' + 'met_percent_wce' + '\t' +
        'exp_chip' + '\t' + 'met_chip' + '\t' + 'non_met_chip' + '\t' + 'met_percent_chip' + '\t' +
        'met_percent_chip - met_percent_wce' + '\n')
    chip_methyl_vs_exp_dict={}
    wce_methyl_vs_exp_dict = {}
    for file in files:
        print(file.name)
        dictionary={}
        for line in file:
            L=line.strip().split()
            gene_name=L[0]
            exp=L[1]
            met=L[2]
            non_met=L[3]
            met_percent=L[4]
            dictionary[gene_name]= (exp,met,non_met,met_percent)

        if 'WCE' not in file.name:
            print(dictionary)
            chip_methyl_vs_exp_dict=dictionary
            print(chip_methyl_vs_exp_dict)
        elif 'WCE' in file.name:
            wce_methyl_vs_exp_dict=dictionary
    for gene in chip_methyl_vs_exp_dict:
        if gene in wce_methyl_vs_exp_dict:
            exp_wce,met_wce,non_met_wce,met_percent_wce=wce_methyl_vs_exp_dict[gene]
            exp_chip,met_chip,non_met_chip,met_percent_chip=chip_methyl_vs_exp_dict[gene]
            chip_nimus_wce_file.write(gene+'\t'+exp_wce+'\t'+met_wce+'\t'+non_met_wce+'\t'+met_percent_wce+'\t'+
                                      exp_chip+'\t'+ met_chip+'\t'+ non_met_chip+'\t'+ met_percent_chip+'\t'+
                                      str(float(met_percent_chip)-float(met_percent_wce))+'\n')

# print(chip_nimus_wce())
def create_plots_expression_vs_methylation():
    chip_methyl_vs_exp_file = open('genes_peeks_4d_RA_H3K36me3_Intersect.ChIP.peaks.CpG.metVSexp.bed')
    wce_methyl_vs_exp_file = open('genes_peeks_4d_RA_H3K36me3_Intersect.WCE.ChIP.peaks.CpG.metVSexp.bed')
    chip_minus_wce_file=open('genes_peeks_4d_RA_H3K36me3_Intersect.WCE.ChIP.peaks.CpG.metVSexp_delta.bed')
    files=[chip_methyl_vs_exp_file,wce_methyl_vs_exp_file,chip_minus_wce_file]
    for file in files:
        exp_y=[]
        met_x=[]
        for line in file.readlines():
            L=line.strip().split()
            if len(L)>=7:
                exp=math.log2(float(L[1])+1)
                # exp=float(L[1])
                exp_y.append(float(exp))
                met_x.append(float(L[9]))
            else:
                exp = math.log2(float(L[1]) + 1)
                # exp = float(L[1])
                exp_y.append(float(exp))
                met_x.append(float(L[4]))
        #plt.hist(met_x,bins=5)
        #plt.show()
        plt.title(file.name)
        plt.plot(met_x,exp_y, 'k.')
        plt.show()

create_plots_expression_vs_methylation()


# print(create_plots(0,1,new_list))


# file=open('genes_peeks_4d_RA_H3K36me3_Intersect.WCE.ChIP.peaks.CpG.metVSexp_delta.bed')
# x=[]
# for line in file.readlines():
#     L=line.strip().split()
#     x.append(L[9])
# print(x)
#
# fig, axs = plt.subplots(1, 2,sharey=True, tight_layout=True)
# axs[0].hist(x)
# axs[1].hist(x, bins=20)
# plt.show()