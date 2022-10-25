#!/usr/bin/env python
# coding: utf-8

# Difexpression tool

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest


# Checking of confidence intervals intersecting
def check_intervals_intersect(first_ci, second_ci):   
    if first_ci[0] < second_ci[0]:
        return second_ci[0] <= first_ci[1] # True or False
    elif first_ci[0] > second_ci[0]:
        return first_ci[0] <= second_ci[1] # True or False
    elif first_ci[0] == second_ci[0]:
        return True
    
# Computing and checking CIs for two tables with genes   
def check_dge_with_ci(first_table, second_table, gene, ci_test_results):
    # dge - differential gene expression
    ft_ci = st.t.interval(confidence=0.95,
    df = len(first_table[gene]) - 1,
    loc = np.mean(first_table[gene]),
    scale = st.sem(first_table[gene])) 
    # NK клетки
    sc_ci = st.t.interval(confidence=0.95,
    df = len(second_table[gene]) - 1,
    loc = np.mean(second_table[gene]),
    scale = st.sem(second_table[gene]))
    ci_test_results.append(check_intervals_intersect(ft_ci, sc_ci))


# Computing z-score
def check_dge_with_ztest(first_table, second_table, gene, z_test_stat, z_test_pvalue):
    # dge - differential gene expression
    # ztest(first_table[gene], second_table[gene])
    #z_test_results.append(ztest(first_table[gene], second_table[gene]))
    z_test_stat.append(np.ravel(ztest(first_table[gene], second_table[gene]))[1]<0.05)
    z_test_pvalue.append(np.round(np.ravel(ztest(first_table[gene], second_table[gene]))[1], 3))

# Computing difference between means of each cell-type

def diff_mean(first_table, second_table, gene, diff_means):
    first_mu = round(np.mean(first_table[gene]), 3)
    second_mu = round(np.mean(second_table[gene]), 3)
    diff_means.append(second_mu - first_mu)


def difexpression_tool(first_cell_type_expressions_path, second_cell_type_expressions_path, save_results_table):
    first_table = pd.read_csv(first_cell_type_expressions_path)
    second_table = pd.read_csv(second_cell_type_expressions_path)
    
    genes1 = np.array(list(first_table.columns), dtype='str')[:-1]
    genes2 = np.array(list(second_table.columns), dtype='str')[:-1]
    genes = []
    ci_test_results = []
    z_test_stat = []
    z_test_pvalue = []
    diff_means = []
    for gene in genes1:
        if gene in genes2:
            
# First column in the final table - gene's name
            genes.append(gene)
    
    
# Second column in the final table - is CIs intersect each other
            check_dge_with_ci(first_table, second_table, gene, ci_test_results)

    
# Third and fourth column in the final table - z-score and p.value
            check_dge_with_ztest(first_table, second_table, gene, z_test_stat, z_test_pvalue)
            #z_test_results = np.matrix(check_dge_with_ztest(first_table, second_table, gene))
            #z_test_stat = np.ravel(z_test_results[:,1])<0.05
            #z_test_pvalue = np.round(np.ravel(z_test_results[:,1]), 3)


# Fifth column in the final table
            diff_mean(first_table, second_table, gene, diff_means)

# Writing output
    results = {
        "gene": genes,
        "ci_test_results": ci_test_results,
        "z_test_results": z_test_stat,
        "z_test_p_values": z_test_pvalue,
        "mean_diff": diff_means
    }
    results = pd.DataFrame(results)
    results[1:].to_csv(save_results_table)
