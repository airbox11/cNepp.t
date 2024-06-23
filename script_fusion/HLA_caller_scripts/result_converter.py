
# coding: utf-8

# ### Convert HLA-typing output files to standard format


# ### Phlat

def convert_phlat(phlat_file):
    import re
    import pandas as pd
    # read Phlat output file in dataframe
    lines = [line.rstrip('\n') for line in open(phlat_file)]    
    lines = [l.split("\t") for l in lines]
    df = pd.DataFrame(lines)
    # set entries in first row as column names and first column as index
    df.columns = df.iloc[0]
    df = df.drop([0])
    df = df.set_index(['Locus'])
    # extract 4 digits from HLA allele name
    matcher = lambda i: re.search(r'.*\*(\d+:\d+):*.*', i).group(1)
    allele1 = df['Allele1'].apply(matcher)
    allele2 = df['Allele2'].apply(matcher)
    # create dictionary of HLA alleles
    d = {'A_1': allele1.at['HLA_A'], 'A_2':allele2.at['HLA_A'], 'B_1':allele1.at['HLA_B'], 'B_2':allele2.at['HLA_B'], 'C_1':allele1.at['HLA_C'], 'C_2':allele2.at['HLA_C'], 'DQA1_1':allele1.at['HLA_DQA1'], 'DQA1_2':allele2.at['HLA_DQA1'], 'DQB1_1':allele1.at['HLA_DQB1'], 'DQB1_2':allele2.at['HLA_DQB1'], 'DRB1_1':allele1.at['HLA_DRB1'], 'DRB1_2':allele2.at['HLA_DRB1']}
    return d


# ### Polysolver

def convert_poly(poly_file):
    import re
    import pandas as pd
    # read Polysolver output file in dataframe
    lines = [line.rstrip('\n') for line in open(poly_file)]    
    lines = [l.split("\t") for l in lines]
    df = pd.DataFrame(lines)
    df = df.set_index([0])
    # extract 4 digits from HLA allele name
    matcher = lambda i: re.search(r'\D+_\D+_(\d+_\d+)_*.*', i).group(1)
    allele1 = df[1].apply(matcher)
    allele2 = df[2].apply(matcher)
    # change format of alleles
    allele1 = allele1.str.replace('_', ':')
    allele2 = allele2.str.replace('_', ':')
    # create dictionary of HLA alleles
    d = {'A_1': allele1.at['HLA-A'], 'A_2': allele2.at['HLA-A'], 'B_1': allele1.at['HLA-B'], 'B_2': allele2.at['HLA-B'], 'C_1': allele1.at['HLA-C'], 'C_2': allele2.at['HLA-C']}
    return d


# ### Optitype

def convert_opti(opti_file):
    import re
    import pandas as pd
    # read Optitype output file in dataframe
    lines = [line.rstrip('\n') for line in open(opti_file)]
    lines = [l.split("\t") for l in lines]
    df = pd.DataFrame(lines).T
    df = df.drop([0,7,8])
    df = df.set_index([0])
    # extract 4 digits from HLA allele name
    matcher = lambda i: re.search(r'.*\*(\d+:\d+):*.*', i).group(1)
    alleles = df[1].apply(matcher)
    # create dictionary of HLA alleles
    d = {'A_1': alleles.loc['A1'], 'A_2': alleles.loc['A2'], 'B_1': alleles.loc['B1'], 'B_2': alleles.loc['B2'], 'C_1': alleles.loc['C1'], 'C_2': alleles.loc['C2']}
    return d

# ### HLA-VBSeq

# HLA-VBSeq gives unfiltered output.
# hlavb_count_occ drops all lines in file with coverage below threshold and returns the filtered output and variables containing the occurence of each hla_locus (A, B, C,...) in the file for further filtering
def hlavb_count_occ(line_list):
    import re
    #filtered_lines = []
    d ={'ind_A':[],'ind_B':[], 'ind_C':[],'ind_DQA1':[],'ind_DQB1':[],'ind_DRB1':[]}
    index = 0
    for line in line_list:
        index+=1
        # filter for coverage and add line index for every allele that passes threshold to list
        if float(line[1]) > 3:
            #filtered_lines.append(line)
            index_real = index -1
            matcher = lambda i: re.search(r'(.*)\*\d+:\d+:*.*', i).group(1)
            al = matcher(line[0])
            d['ind_{}'.format(al)].append(index_real)
    return d
    #return d['ind_A'],d['ind_B'],d['ind_C'],d['ind_DQA1'],d['ind_DQB1'],d['ind_DRB1']


def hlavb_filter(file_path):
    import numpy as np
    import re
    # read HLAVBSeq output file in dataframe
    lines = [line.rstrip('\n') for line in open(file_path)]
    lines = [l.split("\t") for l in lines]
    # get line indexes of every allele that passes threshold
    index_d = hlavb_count_occ(lines)
    #initialize dictionary
    results_d= {'A_1':'', 'A_2':'', 'B_1':'', 'B_2':'', 'C_1':'', 'C_2':'', 'DQA1_1':'', 'DQA1_2':'', 'DQB1_1':'', 'DQB1_2':'', 'DRB1_1':'', 'DRB1_2':''}
    # get allele 1 and 2 for every gene
    for index in index_d:
        # check wich gene (A,B,C,...)
        matcher1 = lambda i: re.search(r'.*_(.*)', i).group(1)
        al = matcher1(index)
        # get 4 digit code
        matcher2 = lambda i: re.search(r'.*\*(\d+:\d+):*.*', i).group(1)
        # no alleles passed threshold
        if len(index_d[index]) == 0:
            results_d['{}_1'.format(al)] = np.nan
            results_d['{}_2'.format(al)] = np.nan
        # 1 allele passed threshold
        elif len(index_d[index]) == 1:
            # homozygous if coverage sufficient
            if float(lines[index_d[index][0]][1]) > 60:
                results_d['{}_1'.format(al)] = matcher2(lines[index_d[index][0]][0])
                results_d['{}_2'.format(al)] = matcher2(lines[index_d[index][0]][0])
            # heterozygous with 1 unknown allele if coverage insufficient
            else:
                results_d['{}_1'.format(al)] = matcher2(lines[index_d[index][0]][0])
                results_d['{}_2'.format(al)] = np.nan
        # 2+ alleles passed threshold
        elif len(index_d[index]) > 1:
            # check if top allele has more than twice the coverage than second candidate --> homozygous
            if float(lines[index_d[index][0]][1]) > 2*float(lines[index_d[index][1]][1]):
                results_d['{}_1'.format(al)] = matcher2(lines[index_d[index][0]][0])
                results_d['{}_2'.format(al)] = matcher2(lines[index_d[index][0]][0])
            # otherwise heterozygous
            else:
                results_d['{}_1'.format(al)] = matcher2(lines[index_d[index][0]][0])
                results_d['{}_2'.format(al)] = matcher2(lines[index_d[index][1]][0])
    return results_d


# ### xHLA

def convert_xhla(xhla_file):
    import re
    # read in xHLA output file
    lines = [line.rstrip('\n') for line in open(xhla_file)]
    lines = [l.split("\"") for l in lines]
    liste = []
    # for every element in line
    for line in lines:
        for l in line:
            if re.match(r'\w+\*\d+:\d+.*', l) :
                liste.append(re.search(r'.*\*(\d+:\d+):*.*', l).group(1))
    #d = {'A_1': liste[0], 'A_2': liste[1], 'B_1': liste[2], 'B_2': liste[3], 'C_1': liste[4], 'C_2': liste[5], 'DPB1_1': liste[6] , 'DPB1_2': liste[7], 'DQB1_1': liste[8], 'DQB1_2': liste[9], 'DRB1_1': liste[10], 'DRB1_2': liste[11]}
    # without DPB1
    d = {'A_1': liste[0], 'A_2': liste[1], 'B_1': liste[2], 'B_2': liste[3], 'C_1': liste[4], 'C_2': liste[5], 'DQB1_1': liste[8], 'DQB1_2': liste[9], 'DRB1_1': liste[10], 'DRB1_2': liste[11]}
    return d


# ### Kourami

def convert_kourami(kourami_file):
    import re
    lines = [line.rstrip('\n') for line in open(kourami_file)]
    lines = [l.split("\t") for l in lines]
    liste = []
    for line in lines:
        for l in line:
            if re.match(r'\w+\*\d+:\d+.*', l) :
                liste.append(re.search(r'.*\*(\d+:\d+):*.*', l).group(1))
    if len(liste) < 12 :
        string = open(kourami_file, 'r').read()
        #print(string)
        if re.match(r'.*DQA1\*\d+:\d+.*', string, re.DOTALL) and re.match(r'.*DQB1\*\d+:\d+.*', string, re.DOTALL):
            d = {'A_1': liste[0], 'A_2': liste[1], 'B_1': liste[2], 'B_2': liste[3], 'C_1': liste[4], 'C_2': liste[5], 'DQA1_1': liste[6] , 'DQA1_2': liste[7], 'DQB1_1': liste[8], 'DQB1_2': liste[9]}
        elif re.match(r'.*DQA1\*\d+:\d+.*', string, re.DOTALL) and re.match(r'.*DRB1\*\d+:\d+.*', string, re.DOTALL):
            d = {'A_1': liste[0], 'A_2': liste[1], 'B_1': liste[2], 'B_2': liste[3], 'C_1': liste[4], 'C_2': liste[5], 'DQA1_1': liste[6] , 'DQA1_2': liste[7], 'DRB1_1': liste[8], 'DRB1_2': liste[9]}
        elif re.match(r'.*DQB1\*\d+:\d+.*', string, re.DOTALL) and re.match(r'.*DRB1\*\d+:\d+.*', string, re.DOTALL):
            d = {'A_1': liste[0], 'A_2': liste[1], 'B_1': liste[2], 'B_2': liste[3], 'C_1': liste[4], 'C_2': liste[5], 'DQB1_1': liste[6] , 'DQB1_2': liste[7], 'DRB1_1': liste[8], 'DRB1_2': liste[9]}
    else :
        d = {'A_1': liste[0], 'A_2': liste[1], 'B_1': liste[2], 'B_2': liste[3], 'C_1': liste[4], 'C_2': liste[5], 'DQA1_1': liste[6] , 'DQA1_2': liste[7], 'DQB1_1': liste[8], 'DQB1_2': liste[9], 'DRB1_1': liste[10], 'DRB1_2': liste[11]}
    
    return d


# ### HLA*PRG:LA

def convert_hlaprg(hlaprg_file):
    import re
    lines = [line.rstrip('\n') for line in open(hlaprg_file)]
    lines = [l.split("\t") for l in lines]
    liste = []
    for line in lines:
        for l in line:
            if re.match(r'\w+\*\d+:\d+.*', l) :
                liste.append(re.search(r'.*\*(\d+:\d+):*.*', l).group(1))
    d = {'A_1': liste[0], 'A_2': liste[1], 'B_1': liste[2], 'B_2': liste[3], 'C_1': liste[4], 'C_2': liste[5], 'DQA1_1': liste[6] , 'DQA1_2': liste[7], 'DQB1_1': liste[8], 'DQB1_2': liste[9], 'DRB1_1': liste[10], 'DRB1_2': liste[11]}
    return d


# ### Create complete resultstable

def combine_results(out_dir):
    import glob
    import pandas as pd
    phlat_path = out_dir + '/PHLAT/*_PHLAT_HLA.sum'
    poly_path = out_dir + '/Polysolver/winners.hla.txt'
    opti_path = out_dir + '/Optitype/20*/20*_result.tsv'
    opti_rna_path = out_dir + '/Optitype_RNA/20*/20*_result.tsv'
    hlavb_path = out_dir + '/HLAVBSeq/*_HLAVBSeq_result_parsed.txt'
    xhla_path = out_dir + '/xHLA/hla-*/*.json'
    kourami_path = out_dir + '/Kourami/*.result'
    hlaprg_path = out_dir + '/HLA_PRG/*/hla/R1_bestguess_G.txt'
    list1 = []
    list2 = []
    if glob.glob(phlat_path):
        list1.append(convert_phlat(glob.glob(phlat_path)[0]))
        list2.append('Phlat') 
    if glob.glob(hlavb_path):
        list1.append(hlavb_filter(glob.glob(hlavb_path)[0]))
        list2.append('HLA-VBSeq')
    if glob.glob(opti_path):
        list1.append(convert_opti(glob.glob(opti_path)[0]))
        list2.append('Optitype') 
    if glob.glob(poly_path):
        list1.append(convert_poly(glob.glob(poly_path)[0]))
        list2.append('Polysolver') 
    if glob.glob(xhla_path):
        list1.append(convert_xhla(glob.glob(xhla_path)[0]))
        list2.append('xHLA') 
    if glob.glob(kourami_path):
        list1.append(convert_kourami(glob.glob(kourami_path)[0]))
        list2.append('Kourami') 
    if glob.glob(hlaprg_path):
        list1.append(convert_hlaprg(glob.glob(hlaprg_path)[0]))
        list2.append('HLA*PRG:LA')
    if glob.glob(opti_rna_path):
        list1.append(convert_opti(glob.glob(opti_rna_path)[0]))
        list2.append('Optitype_RNA')
    res_df = pd.DataFrame(list1)
    res_df.index = list2
    return res_df

def find_consesus(final_df):
    alleles = []
    for allele in ('A', 'B', 'C', 'DQB1', 'DRB1'):
        df_merge = final_df[ allele + '_1'].append(final_df[ allele + '_2'])
        counts = df_merge.value_counts()
        allele1 = counts.keys()[0]
        if len(counts) ==1:
            allele2 = allele1
        elif counts.iloc[0] > 2* counts.iloc[1] :
            allele2 = allele1
        else:
            allele2 = counts.keys()[1]
        alleles.append(allele1)
        alleles.append(allele2)
    return alleles     

def evaluate_performance(final_df, list_alleles):
    import pandas as pd
    #list_alleles = find_consesus(final_df)
    d = {}
    for i in final_df.index:
        count = 0
        correct = 0
        false = 0
        for a in ('A', 'B', 'C'):
            if pd.isnull(final_df.loc[i][count]):
                if pd.isnull(final_df.loc[i][count+1]): 
                    pass
                elif final_df.loc[i][count+1] == list_alleles[count] or list_alleles[count+1] :
                    correct += 1
                else:
                    false += 1
            elif final_df.loc[i][count] == list_alleles[count] :
                correct += 1
                if final_df.loc[i][count+1] == list_alleles[count+1] :
                    correct += 1
                elif pd.isnull(final_df.loc[i][count+1]): 
                    pass
                else:
                    false += 1
            elif final_df.loc[i][count] == list_alleles[count+1] :
                correct += 1
                if final_df.loc[i][count+1] == list_alleles[count] :
                    correct += 1
                elif pd.isnull(final_df.loc[i][count+1]): 
                    pass
                else:
                    false += 1
            elif final_df.loc[i][count] != list_alleles[count] or list_alleles[count+1] :
                false += 1
                if final_df.loc[i][count+1] == list_alleles[count] or list_alleles[count+1] :
                    correct += 1
                elif pd.isnull(final_df.loc[i][count+1]): 
                    pass
                else:
                    false += 1
            count += 2
        d.update({i+'_correct' : correct, i+'_false' : false})
    return d
 
def input_netMHCI(final_df):
    input_netMHCI = ''
    for allele in ('A', 'B', 'C'):
        df_merge = final_df[ allele + '_1'].append(final_df[ allele + '_2'])
        counts = df_merge.value_counts()
        allele1 = counts.keys()[0]
        if len(counts) ==1:
            allele2 = allele1
        elif counts.iloc[0] >  2* counts.iloc[1] :
            allele2 = allele1
        else:
            allele2 = counts.keys()[1]
        if len(input_netMHCI) > 0:
            input_netMHCI = input_netMHCI + ',HLA-' + allele + allele1 + ',HLA-' + allele + allele2
        else:
            input_netMHCI = input_netMHCI + 'HLA-' + allele + allele1 + ',HLA-' + allele + allele2
    return input_netMHCI

def input_netMHCII(final_df):
    input_netMHCII = ''
    for allele in ('DQB1', 'DRB1'):
        df_merge = final_df[ allele + '_1'].append(final_df[ allele + '_2'])
        counts = df_merge.value_counts()
        allele1 = counts.keys()[0].replace(':','')
        if len(counts) ==1:
            allele2 = allele1
        elif counts.iloc[0] >  2* counts.iloc[1] :
            allele2 = allele1
        else:
            allele2 = counts.keys()[1].replace(':','')
        if len(input_netMHCII) > 0:
            input_netMHCII = input_netMHCII + ',' + allele + '_' + allele1 + ',' + allele + '_' + allele2
        else:
            input_netMHCII = allele + '_' + allele1 + ',' + allele + '_' + allele2
    return input_netMHCII      
    
def input_netMHC(final_df, out_dir):
    import glob
    with open(glob.glob(out_dir)[0] + '/netMHCI_input.txt', "w") as text_file:
        text_file.write(input_netMHCI(final_df))
    with open(glob.glob(out_dir)[0] + '/netMHCII_input.txt', "w") as text_file:
        text_file.write(input_netMHCII(final_df))

def results_with_consensus(final_df, out_dir):
    import glob
    import pandas as pd
    d = {}
    for allele in ('A', 'B', 'C', 'DQB1', 'DRB1'):
        df_merge = final_df[ allele + '_1'].append(final_df[ allele + '_2'])
        counts = df_merge.value_counts()
        allele1 = counts.keys()[0]
        d.update({ allele + '_1' : allele1 })
        if len(counts) ==1:
            allele2 = allele1
            d.update({ allele + '_2' : allele2 })
        elif counts.iloc[0] >  2* counts.iloc[1] :
            allele2 = allele1
            d.update({ allele + '_2' : allele2 })
        else:
            allele2 = counts.keys()[1]
            d.update({ allele + '_2' : allele2 })
    cons = pd.DataFrame(d, index=['Consensus',])
    out_df = pd.concat([final_df, cons])
    out_csv = out_df.to_csv()
    with open(glob.glob(out_dir)[0] + '/HLA_results_with_consensus.csv', "w") as text_file:
        text_file.write(out_csv)
        