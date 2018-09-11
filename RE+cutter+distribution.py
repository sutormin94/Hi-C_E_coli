
###############################################
###Based on Monica Guo script####
##Dmitry Sutormin, 2018##
##E. coli Hi-C project##

#Script analysis the distribution of restriction sites across the genome.
#It takes genome FASTA as input and a dictionary of restrictases {name: site sequence}.
#Plots distribution of restriction fragments lengths and distribution of the restriction sites numbers fall into some bins.
###############################################


#######
#Packages to be imported.
#######

import os
import regex
from Bio.Seq import Seq
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors

#######
#Variables to be defined.
#######

Genome_c_cres_path="C:\Sutor\science\Caulobacter_Topo-Seq\Genome\C_crescentus_NC_011916.1.fasta"
Genome_e_coli_path="C:\Sutor\science\DNA-gyrase\Genomes\E_coli_w3110_G_Mu.fasta"

Path_out="C:\Sutor\science\E_coli_Hi-C\RE_to_use\RE_analysis\\"

# 'BglII': 'AGATCT', bad has 6617.54/5034.0 cuts. use as pos ctrl
# 'AseI': 'AT|TA|AT' 2 bp 5' overhang. Not methylation sensitive. 2489.78/1354.0 cuts. R0526M = 50 kU Great choice
# 'BssHII': 'GCGCGC', Bad, digest at 50C. Can buy PauI or PteI from Lifetech, but may not be ideal. PauI is best - 4 bp 5' overhang, not methyl sensitive, 1871.85/1143.0 cuts. G|CGCG|C 
# 'DpnII': '|GATC|', methylation dependent, could use a mix of DpnI and DpnII? 
# 'HpaII': 'C|CG|G', 2 bp 5' overhang. Not methylation sensitive, Inhibited by salt concentrations > 50 mM KCl. has 50 kU size 190.91/122.0

REs = {'BglII': 'AGATCT', 'AseI': 'ATTAAT', 'PauI': 'GCGCGC', 'NcoI': 'CCATGG', 'EcoRI': 'GAATTC', 'DpnI/II': 'GATC', 'HpaII': 'CCGG', 'HhaI': 'GCGC'}


#######
#FASTA parser.
#######

def obtain_seq(seq_path):
    seq_oi=open(seq_path, 'r')
    for record in SeqIO.parse(seq_oi, "fasta"):
        sequence=str(record.seq)
    seq_oi.close()      
    return sequence, record.id

#######
#Looks for restriction sites in the genome.
#######

def calculateREs(RE, REsite, genome):
    #Finds all sites
    allSites = regex.finditer(REsite, genome, overlapped=True)
    sites = [] #List with sites coordinates
    for i in allSites:
        sites.append(i.start())
    #List with restriction fragments length
    binsize = []
    for i in range(len(sites)-1):
        binsize.append(sites[i+1]-sites[i])
    re_fragment_mean=np.around(np.mean(binsize))
    re_fragment_median=np.around(np.median(binsize))
    print("{} mean {}, median {}".format(RE, re_fragment_mean, re_fragment_median))
    #Plot distribution of restriction fragments length
    #plt.figure()
    #plt.hist(binsize)
    #plt.title('{} mean {}, median {}'.format(RE, np.around(np.mean(binsize)), np.around(np.median(binsize))))
    #plt.xlabel('size of {} fragment'.format(RE))
    #plt.close()
    return sites, binsize, re_fragment_mean, re_fragment_median

#######
#Finds bins (long genome regions) which are not covered with restriction sites
#######

def combine_two_RE(Name_ar, RE_data_dict):
    all_sites=[]
    name_comb=''
    for name in Name_ar:
        all_sites+=RE_data_dict[name][0]
        name_comb+=name+'/'
    name_comb=name_comb.rstrip('/')
    all_sites=sorted(all_sites)
    #List with restriction fragments length
    binsize = []
    for i in range(len(all_sites)-1):
        binsize.append(all_sites[i+1]-all_sites[i])
    re_fragment_mean=np.around(np.mean(binsize))
    re_fragment_median=np.around(np.median(binsize))
    print("{} mean {}, median {}".format(name_comb, re_fragment_mean, re_fragment_median))  
    return all_sites, binsize, re_fragment_mean, re_fragment_median

#######
#Finds bins (long genome regions) which are not covered with restriction sites
#######

def calcEmptys(RE, REsites_list, genome, binsize = 10000):
    #Calculates numbers of restriction sites falling into bins (evenly distributed genome regions).
    bins = range(0, len(genome), binsize)
    Sites_in_bins = []
    for i in range(len(bins)):
        if i == len(bins)-1:
            Sites_in_bins.append(sum([bins[i] <= x < len(genome) for x in REsites_list]))
        else:
            Sites_in_bins.append(sum([bins[i] <= x < bins[i+1] for x in REsites_list]))
    print("Number of {} sites {}".format(RE, len(REsites_list)))
    print("Sum of {} sites in bins {}".format(RE, sum(Sites_in_bins)))
    #Finds bins that do not contain restriction sites
    Bins_no_site=np.where(Sites_in_bins <= np.array(Sites_in_bins).min())[0] #Get empty bens indexes
    Bins_no_site_ratio=len(Bins_no_site)/float(len(Sites_in_bins))
    print("Number of bins without {} sites {}".format(RE, len(Bins_no_site)))
    print("Ratio of bins without {} sites {}\n".format(RE, Bins_no_site_ratio))
    #Plot distribution of the number of cuts that fall into bins
    #plt.figure()
    #plt.hist(Sites_in_bins)
    #plt.title(str(RE))
    #plt.xlabel('{} cuts per bin'.format(RE))
    #plt.close()
    return Sites_in_bins, Bins_no_site, len(Bins_no_site), Bins_no_site_ratio

#######
#Sorts restrictases by their performance: num of sites, num of bins empty, ratio of empty bins
#######

def rank_RE(RE_data):
    #Rank according to individual parameters
    RE_sorted_by_fragment_mean=sorted(RE_data, key=lambda tup: tup[3]) #Lower length - better
    RE_sorted_by_fragment_median=sorted(RE_data, key=lambda tup: tup[4]) #Lower length - better
    RE_sorted_by_num_sites=sorted(RE_data, key=lambda tup: tup[5], reverse=True) #More sites - better
    RE_sorted_by_num_bins_empty=sorted(RE_data, key=lambda tup: tup[8]) #Less empty - better
    RE_sorted_by_num_bins_empty_ratio=sorted(RE_data, key=lambda tup: tup[9]) #Less ratio - better 
    #Combining the parameters
    RE_ranked=[]
    for RE_info in RE_sorted_by_num_sites:
        mean_index=RE_sorted_by_fragment_mean.index(RE_info)
        median_index=RE_sorted_by_fragment_median.index(RE_info)
        ns_index=RE_sorted_by_num_sites.index(RE_info)
        be_index=RE_sorted_by_num_bins_empty.index(RE_info)
        ber_index=RE_sorted_by_num_bins_empty_ratio.index(RE_info)
        RE_ranked.append([RE_info[0], RE_info[3], RE_info[4], RE_info[5], RE_info[8], RE_info[9]]+[mean_index, median_index, ns_index, be_index, ber_index, mean_index+median_index+ns_index+be_index+ber_index])
    RE_ranked=sorted(RE_ranked, key=lambda tup: tup[11])
    return RE_ranked

#######
#Plotting the distribution of RE fragments length
#######

def plot_fragments_length(RE_ranked, RE_data_dict, path_out):
    #Set number of plots
    Num_plots=len(RE_ranked)
    
    #Plot distribution of restriction fragments length
    if Num_plots%4==0:
        fig1, plots=plt.subplots(Num_plots//4,4,figsize=(15,15), dpi=100)
    else:
        fig1, plots=plt.subplots((Num_plots//4)+1,4,figsize=(15,15), dpi=100)
    i=0
    iter_plot=plots.reshape(-1)
    while i<Num_plots:
        RE_name=RE_ranked[i][0]
        fragsize=RE_data_dict[RE_name][1]
        iter_plot[i].hist(fragsize, edgecolor='black', linewidth=1)
        iter_plot[i].set_title(RE_name, size=22)
        iter_plot[i].set_xlabel('Size of fragments', size=18)  
        iter_plot[i].set_ylabel('Number of fragments, bp', size=18)
        iter_plot[i].tick_params(axis='both', labelsize=14)
        iter_plot[i].annotate('mean {} bp\nmedian {} bp'.format(int(np.mean(fragsize)), int(np.median(fragsize))), xytext=(0.2, 0.8), textcoords='axes fraction', xy=(0, 0), size=16.5)
        i+=1
    
    #Remove empty subplots
    for elem in iter_plot[Num_plots:]:
        plt.delaxes(elem)
    
    plt.tight_layout()
    plt.savefig(path_out+"E_coli_re_analysis_fragment_length.png", figsize=(15,10), dpi=400)
    plt.close()
    return

#######
#Plotting the distribution of the number of cuts in genome bins
#######

def plot_cut_in_bin(RE_ranked, RE_data_dict, path_out):
    #Set number of plots
    Num_plots=len(RE_ranked)
    
    #Plot distribution of restriction fragments length
    if Num_plots%4==0:
        fig1, plots=plt.subplots(Num_plots//4,4,figsize=(15,15), dpi=100)
    else:
        fig1, plots=plt.subplots((Num_plots//4)+1,4,figsize=(15,15), dpi=100)
    i=0
    iter_plot=plots.reshape(-1)
    while i<Num_plots:
        RE_name=RE_ranked[i][0]
        catnum=RE_data_dict[RE_name][2]
        iter_plot[i].hist(catnum, edgecolor='black', color='#ff862f', linewidth=1)
        iter_plot[i].set_title(RE_name, size=22)
        iter_plot[i].set_xlabel('Number of cuts', size=18)  
        iter_plot[i].set_ylabel('Number of bins', size=18)
        iter_plot[i].tick_params(axis='both', labelsize=14)
        print('mean '+str(int(np.mean(catnum)))+' cuts\nmedian '+str(int(np.median(catnum)))+' cuts')
        iter_plot[i].annotate('mean '+str(int(np.mean(catnum)))+' cuts\nmedian '+str(int(np.median(catnum)))+' cuts', xytext=(0.4, 0.8), textcoords='axes fraction', xy=(0, 0), size=16.5)
        i+=1
    
    #Remove empty subplots
    for elem in iter_plot[Num_plots:]:
        plt.delaxes(elem)
    
    plt.tight_layout()
    plt.savefig(path_out+"E_coli_re_analysis_cut_in_bin.png", figsize=(15,10), dpi=400)
    plt.close()
    return

#######
#Plotting the number of empty bins or ratio of empty bins
#######

def empty_num_and_ratio_plot(RE_ranked, RE_data_dict, path_out):
    position=[]
    names=[]
    value_num=[]
    value_ratio=[]
    for i in range(len(RE_ranked)):
        position.append(i)
        names.append(RE_ranked[i][0])
        value_num.append(RE_ranked[i][4])
        value_ratio.append(RE_ranked[i][5])
    fig2, plots=plt.subplots(2,1,figsize=(12,12), dpi=100)
    #Number of empty bins
    plots[0].bar(position, value_num, 0.9, color='#6797ff', linewidth=1, edgecolor='black', align='edge')
    plots[0].set_title('Restrictases and the number of empty bins', size=22)
    plots[0].set_ylabel('Number of empty 10000 bp bins', size=18)
    plots[0].set_ylim(0, max(value_num)+10)
    plots[0].set_xticks(position, minor=False)
    plots[0].set_xticklabels(names)
    plt.setp(plots[0].set_xticklabels(names), rotation=45, fontsize=18)  
    plots[0].yaxis.grid(True, linewidth=1, linestyle='--', color='black')
    plots[0].tick_params(axis='both', labelsize=14)
    for i in range(len(value_num)):
        plots[0].annotate(str(value_num[i]), xytext=(i+0.2, value_num[i]+4), textcoords='data', xy=(0, 0), size=13, weight="bold")
    #Ratio of empty bins
    plots[1].bar(position, value_ratio, 0.9, color='#7eff67', linewidth=1, edgecolor='black', align='edge')
    plots[1].set_title('Restrictases and the ratio of empty bins', size=22) 
    plots[1].set_ylabel('Ratio of empty 10000 bp bins', size=18)
    plots[1].set_ylim(0, max(value_ratio)+0.05)
    plots[1].set_xticks(position, minor=False)
    plots[1].set_xticklabels(names)
    plt.setp(plots[1].set_xticklabels(names), rotation=45, fontsize=18)  
    plots[1].yaxis.grid(True, linewidth=1, linestyle='--', color='black')    
    plots[1].tick_params(axis='both', labelsize=14)  
    for i in range(len(value_ratio)):
        plots[1].annotate(str(round(value_ratio[i], 3)), xytext=(i, value_ratio[i]+0.01), textcoords='data', xy=(0, 0), size=13, weight="bold")    
    plt.tight_layout()
    plt.savefig(path_out+"E_coli_re_analysis_empty_bins.png", figsize=(12,12), dpi=400)
    plt.close()
    return

#######
#Wrapps functions for RE pairs combining
#######

def RE_pair_combining(Name_ar, RE_data_dict, genome, binsize, RE_data):
    #New name - combination
    name_comb=''
    for name in Name_ar:
        name_comb+=name+'/'
    name_comb=name_comb.rstrip('/')    
    #Combined sites
    Combined_sites=combine_two_RE(Name_ar, RE_data_dict)
    #Combined bins
    Combined_bins_analysis=calcEmptys(name_comb, Combined_sites[0], genome[0], binsize)
    RE_data.append([name_comb]+list(Combined_sites)+[len(Combined_sites[0])]+list(Combined_bins_analysis))
    RE_data_dict[name_comb]=[Combined_sites[0], Combined_sites[1], Combined_bins_analysis[0]]
    return

#######
#Wrapps all the functions
#######

def wrapper(REs_dict, genome_path, binsize, path_out):
    genome=obtain_seq(genome_path)
    RE_data=[]
    RE_data_dict={}
    for RE, REsite in REs_dict.items():
        RE_sites_list=calculateREs(RE, REsite, genome[0])
        Bins_analysis=calcEmptys(RE, RE_sites_list[0], genome[0], binsize)
        #RE name [0], Sites coordinates [1], RE fragments length [2], RE fragment mean [3], RE fragment median [4], Num of sites [5], 
        #Num sites in bins [6], Bins no sites [7], Num of empty bins [8], Ratio of bins no sites [9]       
        RE_data.append([RE]+list(RE_sites_list)+[len(RE_sites_list[0])]+list(Bins_analysis)) 
        RE_data_dict[RE]=[RE_sites_list[0], RE_sites_list[1], Bins_analysis[0]] #Sites coords, RE fragments length
    
    #Probing different combinations of RE
    RE_pair_combining(['AseI', 'BglII'], RE_data_dict, genome, binsize, RE_data)
    RE_pair_combining(['EcoRI', 'NcoI'], RE_data_dict, genome, binsize, RE_data)
    RE_pair_combining(['PauI', 'EcoRI'], RE_data_dict, genome, binsize, RE_data)
    RE_pair_combining(['NcoI', 'AseI'], RE_data_dict, genome, binsize, RE_data)
    RE_pair_combining(['PauI', 'AseI'], RE_data_dict, genome, binsize, RE_data)
    RE_pair_combining(['PauI', 'AseI', 'BglII'], RE_data_dict, genome, binsize, RE_data)
    RE_pair_combining(['EcoRI', 'NcoI', 'AseI'], RE_data_dict, genome, binsize, RE_data)
    
    #Ranking RE and combinations of RE
    RE_ranked=rank_RE(RE_data)
    print(RE_ranked)
    
    plot_cut_in_bin(RE_ranked, RE_data_dict, path_out)
    plot_fragments_length(RE_ranked, RE_data_dict, path_out)
    empty_num_and_ratio_plot(RE_ranked, RE_data_dict, path_out)
    return

#wrapper(REs, Genome_c_cres_path, 10000)
wrapper(REs, Genome_e_coli_path, 10000, Path_out)

print('Script ended its work succesfully!') 