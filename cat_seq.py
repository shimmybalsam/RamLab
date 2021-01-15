import glob
import numpy
import csv
import os
import statistics
import matplotlib.pyplot as plt
from pathlib import Path
# p = Path('.')
# all_files=list(p.glob('cat_Sequential_RA_0413_0418/*/*.bed'))

p=Path('C:/Users\חוה רוט/Google Drive/.')
all_files=list(p.glob('cat_Sequential_RA_0413_0418/*/*.bed'))
all_fc_files=list(p.glob('cat_Sequential_RA_0413_0418/*.bed'))
print(len(all_files))
print(all_files)
print(len(all_fc_files))
print(all_fc_files)
all_new=Path('.')
all_new_path=list(all_new.glob('new/*.bed'))
print(all_new_path)

for fc_file in all_fc_files:
    print('fc_file.name',fc_file.name)

import re
import math

# fc_files=Path()
# print(fc_files)

# all_new_files=list(pn.open())
# print(all_new_files)

# all_new_files=
# print(all_new_files)
#
# print(all_new_files)

# print(help(csv))
# files_exons=[]
# for f in all_files:
#     if 'exon' in f.name:
#         files_exons.append(f)
# print(files_exons)


# chip=0
# wce=1

def gene_coordinate(file_coordinates):
    file=open(file_coordinates).readlines()
    all_coordinates={}
    for line in file:
        L = line.strip().split()
        chrom=L[5]
        start=L[6]
        stop=L[7]
        name=L[8]
        if name not in all_coordinates:
            all_coordinates[name]=(chrom,start,stop)
    print(all_coordinates)
    return all_coordinates
# print(gene_coordinate())

def fc_peeks_coordinates(fc_file,file_coordinates):
    # for file in all_files:
    #     print()
    #     if file.parent.name in fc_file.name and 'exons' not in file.name and 'WCE' not in file.name:
    #         file_coordinates=file
    #         print('file',file)
    all_coordinates=gene_coordinate(file_coordinates)
    file_fc=open(fc_file)
    fc_peeks = {}
    for line in file_fc.readlines():
        L = line.strip().split()
        chrom,start,stop=all_coordinates[L[0]]
        fc_peeks[L[0]]=(chrom,start,stop)
    return fc_peeks
# print(fc_peeks_coordinates())

def gene_names_coordinates():
    mm9_refseq_file=open('mm9_RefSeq')
    # new_file = open('peeks_coordinates_Intersect.ChIP.peaks.CpG.bed', 'a')
    mm9_dict={}
    for line in mm9_refseq_file.readlines():
        L = line.strip().split()
        chrom=L[2]
        gene_name=L[12]
        txstart=L[4]
        txstop=L[5]
        mm9_dict[gene_name]=(chrom,txstart,txstop)
    return mm9_dict

def fc_peeks_files(all_fc_files):
    for fc_files in all_fc_files:
        file_coordinates=None
        for file in all_files:
            print()
            if file.parent.name in fc_files.name and 'exons' not in file.name and 'WCE' not in file.name:
                file_coordinates = file
                print('file', file)
        fc_peeks = fc_peeks_coordinates(fc_files,file_coordinates)
        fc_file = fc_files.open().readlines()
        # fc_peeks=fc_peeks_coordinates(fc_file)
        new_fc_file= open('genes_fc_'+fc_files.name , 'a')
        mm9_dict=gene_names_coordinates()
        genes_in_peeks={}
        for k, i in fc_peeks.items():
            peek_name=k
            chr=i[0]
            start=int(i[1])
            end = int(i[2])
            counter = 0
            for mm9_k,mm9_i in mm9_dict.items():
                gene=mm9_k
                chrom_mm9=mm9_i[0]
                start_mm9=mm9_i[1]
                stop_mm9 = mm9_i[2]
                original = gene
                if chrom_mm9 == chr:
                    if (int(start_mm9)<=int(start) and int(stop_mm9)>=int(start)) or (int(start_mm9)<=int(end) and int(stop_mm9)>=int(end)):
                        if gene in genes_in_peeks:

                            gene=gene+'.'+str(counter)
                            genes_in_peeks[gene] = (original,peek_name, chr, start, end)
                            counter+=1
                        else:
                            genes_in_peeks[gene] = (gene,peek_name,chr,start,end)
        print(genes_in_peeks)
        for k,i in genes_in_peeks.items():
            new_fc_file.write(str(str(k)+'\t'+str(i[0])+'\t'+str(i[1])+'\t'+str(i[2])+'\t'+str(i[3])+'\t'+str(i[4])+'\t'+'\n'))
fc_peeks_files(all_fc_files)

def remove_duplicate_exons(all_files):
    for files in all_files:
        all=set()
        if 'exon' in files.name:
            path=Path(files.parent)
            print(files.parent)
            name='no_duplicate_'+ files.name
            print(name)
            new=open(path/name,'a')
            file = files.open().readlines()
            lst = []
            # all = {}
            for line in file:
                # print(line)
                L = line.strip().split()
                if L[0:3] not in lst:
                    # print('lst')
                    new.writelines(line)
                    lst.append(L[0:3])
                else:
                    lst.append(L[0:3])
                # print(lst)
                # str='t'
                # print(line.index(str))
                # if line[0:line.index't0')] not in lst:
                # L = line.strip().split()
                # T=tuple(L)
                # all.add(T)
                # for i in all:
                #     new.writelines(line)
                #     print(lst)
                #     lst.append(line)
            # print(all)
        # print(all)

# print(remove_duplicate_exons(all_files))

# new_all_files=[file for file in all_files if 'no_duplicate' in file.name]
# print('new_all_files',new_all_files)


def exons_all_cpg(all_files):
    """ all_files= only 'exon' in the file name from no_duplicate files.
    marge all exons cpg according to exons"""

    for files in all_files:
        print(files.name)
        new=open('new/new_'+files.parent.name+files.name,'a')
        print(new)
        file=files.open().readlines()
        lst = []
        all_cpg = {}
            # all_methyl=0
            # all_not_methyl=0

        for line in file:
            L = line.strip().split()
            # if 'exon' in files.name:
                # transcript_name = L[15]
                # name=transcript_name
            # else:
            #     pick_name=L[8]
            #     name=pick_name
            chrom=L[12]
            start_exon=L[13]
            stop_exon=L[14]
            name=(chrom,start_exon,stop_exon)
            all_methyl = int(L[3])
            all_not_methyl = int(L[4])
            reads = 1
            if name in all_cpg:
                methyl, no_methyl, reads = all_cpg[name]
                all_cpg[name] = (int(methyl) + int(all_methyl), int(no_methyl) + int(all_not_methyl), reads + 1)
            else:
                all_cpg[name] = (all_methyl, all_not_methyl, 1)
        for k,i in all_cpg.items():
            new.write(str(str(k[0])+'\t'+str(k[1])+'\t'+str(k[2])+'\t'+str(i[0])+'\t'+str(i[1])+'\t'+str(i[2])+'\t'+(str(float(i[0]/(i[0]+i[1]))))+'\t'+'\n'))
        lst.append(all_cpg)
        print(all_cpg)
        print()
        new.close()
# print(exons_all_cpg(new_all_files))
print()



def files_all_cpg(all_files):
    """ merge according to picks """
    for files in all_files:
        print(files.name)
        new=open('new/new_'+files.parent.name+files.name,'a')
        print(new)
        file=files.open().readlines()
        lst = []
        all_cpg = {}
            # all_methyl=0
            # all_not_methyl=0

        for line in file:
            L = line.strip().split()
            # if 'exon' in files.name:
                # transcript_name = L[15]
                # name=transcript_name
            # else:
            #     pick_name=L[8]
            #     name=pick_name
            name=L[8]
            all_methyl = int(L[3])
            all_not_methyl = int(L[4])
            reads = 1
            if name in all_cpg:
                methyl, no_methyl, reads = all_cpg[name]
                all_cpg[name] = (int(methyl) + int(all_methyl), int(no_methyl) + int(all_not_methyl), reads + 1)
            else:
                all_cpg[name] = (all_methyl, all_not_methyl, 1)
        for k,i in all_cpg.items():
            new.write(str(str(k)+'\t'+str(i[0])+'\t'+str(i[1])+'\t'+str(i[2])+'\t'+(str(float(i[0]/(i[0]+i[1]))))+'\t'+'\n'))
        lst.append(all_cpg)
        print(all_cpg)
        print()
        new.close()
# print(files_all_cpg(all_files))
print()



def add_pick_length(all_files,all_new_files):
    for files in all_files:
        print(files)
        for new_files in all_new_files:
            print(new_files)
            if files.parent.name+files.name in new_files.name:
                new_file=new_files.open()
                print(new_files.name)
                file = files.open().readlines()
                lst = []
                all_len = {}
                for line in file:
                    L = line.strip().split()
                    # if 'exon' in files.name:
                    # transcript_name = L[15]
                    # name=transcript_name
                    # else:
                    #     pick_name=L[8]
                    #     name=pick_name
                    name = L[8]
                    start=L[6]
                    stop=L[7]
                    length=L[9]
                    if name in all_len:
                        continue
                    else:
                        all_len[name] = (start, stop, length)
                for new_line in new_file.readlines():
                    new_L = new_line.strip().split()
                    print(new_L)
                    name=new_L[0]
                    start,stop,length=all_len[name]
                    new_file[new_line].write('\t'+start+'\t'+stop+'\t'+length+'\n')


# all_add_pick_new_files=[file for file in all_new_path if 'exons' not in file.name]
# print('all_add_pick_new_files',all_add_pick_new_files)
# all_add_pick_files =[file for file in all_files if 'exons' not in file.name]

# print('all_add_pick_files',all_add_pick_files)
#
# for files in all_add_pick_files:
#     print()
#     for new_files in all_add_pick_new_files:
#         if files.parent.name + files.name in new_files.name:
#             print('files.parent.name + files.name',files.parent.name + files.name)
#             print('new_files.name', new_files.name)
# print(add_pick_length(all_add_pick_files[0:2],all_add_pick_new_files[0:2]))

def single_cpg(all_files):
    """ all files without exons"""
    for files in all_files:
        if 'exons' not in files.name:
            print(files.name)
            new=open('new/single_cpg_'+files.parent.name+files.name,'a')
            print(new)
            file=files.open().readlines()
            lst = []
            all_cpg = {}
            for line in file:
                L = line.strip().split()
                if int(L[3])+int(L(4))>=5:
                    new.writelines(line)


# new_list = [file for file in all_new_path if 'H3K36me3' in file.name and 'WCE' not in file.name]
def create_new_list(all_new_path):
    new_list = [file for file in all_new_path if 'H' in file.name and 'WCE' not in file.name]
    print('new_list',new_list)
    for i in new_list:
        if 'exon' in i.name and 'no_duplicate' not in i.name:
            new_list.remove(i)
    print(new_list)
    print(len(new_list))
    return(new_list)
# print(create_new_list(all_new_path))




def all_data_lst(new_list,cpg,exons):
    x = []
    y = []
    for files in new_list[cpg:exons+1]:

        # print(files.name)
        # print(new_list.index(files))
        #
        if new_list.index(files) % 2 == 0: # cpg files
            print(files.name)
            print(new_list.index(files))
            print(new_list.index(files) % 2)
            file = files.open().readlines()
            for line in file:
                L = line.strip().split()
                # print(L)
                x.append(float(L[4]))
                # print('x',x)
        elif new_list.index(files)%2 != 0:  # exons files
            print(files.name)
            print(new_list.index(files))
            print(new_list.index(files) % 2)
            file = files.open().readlines()
            for line in file:
                L = line.strip().split()
                # y.append(float(L[4]))
                y.append(float(L[6]))
    return (x,y)


def create_plots(cpg,exons,new_list):
    if exons<=len(new_list):
        all_data= all_data_lst(new_list,cpg,exons)

        print(statistics.mean(all_data[0]))
        print(statistics.median(all_data[0]))
        print(statistics.mean(all_data[1]))
        print(statistics.median(all_data[1]))
        plt.rc('font', size=8)
        labels = [(new_list[cpg].name[4:10]+'\n'+new_list[cpg].name[10:32]+'\n'+new_list[cpg].name[32:-4])
            ,(new_list[exons].name[4:10]+'\n'+new_list[exons].name[10:32]+'\n'+new_list[exons].name[32:-4])]

        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))
        print()
        print(fig)
        # print(axes[0])
        # print(axes[1])
        #
        # # rectangular box plot
        # bplot1 = axes.boxplot(all_data,
        #                          vert=True,  # vertical box alignment
        #                          patch_artist=True)  # fill with color
        # # labels=labels)  # will be used to label x-ticks
        # print(bplot1)
        # axes.set_title('Rectangular box plot')

        # notch shape box plot
        bplot2 = axes.boxplot(all_data,
                                 notch=True,  # notch shape
                                 vert=True,  # vertical box alignment
                                 patch_artist=True, showmeans=True, meanline=True,  # fill with color
        labels=labels)  # will be used to label x-ticks
        axes.set_title('Notched box plot')

        # fill with colors
        colors = ['c', 'lawngreen']
        # for bplot in (bplot1, bplot2):
        for patch, color in zip(bplot2['boxes'], colors):
            patch.set_facecolor(color)

        # for ax in axes:
        axes.yaxis.grid(True)
        axes.set_xticks([y + 1 for y in range(len(all_data))],)
            # plt.semilogy(t, np.exp(-t / 5.0))
            # ax.set_xlabel('xlabel')
            # ax.set_ylabel('ylabel')


        plt.show()
        # plt.savefig('new_list[cpg].name.jpg')
        # cpg+=2
        # exons+=2
        create_plots(cpg+2,exons+2,new_list)
# print(create_plots(0,1,new_list))

    # print('all',len(all[-1]))
# print(glob.glob("c:\Google Drive\cat_Sequential_RA_0413_0418/*/*.BED"))
# print("https://drive.google.com/drive/folders/1px1XwqKYQ23iVpqPpE6IP3etvxo2wXBV/*.BED")

# https://drive.google.com/drive/folders/1px1XwqKYQ23iVpqPpE6IP3etvxo2wXBV

# import os
# arr = os.open('cat_Sequential_RA_0413_0418')
# for a in arr:
#
#    print(a.path)
# print(arr)