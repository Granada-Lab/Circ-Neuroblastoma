###############################################################################################
# SCRIPT FOR ANALYSIS PIPELINE, COLLABORATION WITH CHARITE-BERLIN, NEUROBLASTOMA # March 2025 #
###############################################################################################

# SysBio Group, DLSM, FSTC, University of Luxembourg
# PhD candidate Jeff Didier
# Dr Sebastien De Landtsheer
# Prof. Dr-Ing Thomas Sauter

# Running time cold: 10 minutes

#########################
# ## LIBRARIES IMPORTS ##

import os
import random
import itertools

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from adjustText import adjust_text
from scipy.stats import spearmanr
from functions.plot_functions import plot_LDA
from functions.data_functions import get_LDA_metrics, LDA_loocv
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

#########################################################################
# ## CREATE FIX VARIABLES, CREATE RESULTS FOLDER, AND COLLECT RAW DATA ##

# Variables
dpi = 300
seed = 42
formats = ['svg', 'png']
res_folder = 'results'

ccle_exp = 'CCLE_expression_from_CCLE_2022_11_22.csv'
circ_data = 'neuroblastoma_circadian_data_v1.xlsx'
gene_lists = 'neuroblastoma_gene_lists.xlsx'
sample_info = 'sample_info.csv'
drug_sensitivities = 'CCLE_NP24.2009_Drug_data_2015.02.24.csv'

# Results
if os.path.isdir(f'./{res_folder}') is False:
    os.mkdir(f'./{res_folder}')

# Raw data
circ_score = pd.read_excel(f"../data/{circ_data}", sheet_name='circ_score')
circ_values = pd.read_excel(f"../data/{circ_data}", sheet_name='circ_values_selection')

circ_genes = pd.read_excel(f"../data/{gene_lists}", sheet_name='clock_genes')
neuro_genes = pd.read_excel(f"../data/{gene_lists}", sheet_name='neuroblastoma_genes')

data_exp = pd.read_csv(f"../data/{ccle_exp}")
drug_sens = pd.read_csv(f"../data/{drug_sensitivities}")
data_info = pd.read_csv(f"../data/{sample_info}").set_index('CCLE_Name')

# Random number generator
random.seed(seed)
np.random.seed(seed)

###############################################################################################
# ## Cleaning CCLE for clock and neuroblastoma genes, drug sensitivity and expression values ##

# Collect cell line names
celllines = circ_score['Cellline']

# 3 celllines need adaptations to find their ACH-code to map with CCLE expression: Lan5, SKNBE, Kelly
celllines = [x.upper() for x in celllines]  # Kelly -> KELLY, Lan5 -> LAN5
celllines[2] = 'SKNBE2'  # SKNBE -> SKNBE2

# Get the DepMapID of each cell line (ACH-code) (only 9 of 11 found)
DepMapId = [data_info[data_info['stripped_cell_line_name'] == x].loc[:, 'DepMap_ID'].to_numpy()[0]
            if len(data_info[data_info['stripped_cell_line_name'] == x].loc[:, 'DepMap_ID'].to_numpy()) != 0 else
            '' for x in celllines]

# Get the corresponding expression values (only 8 of 9 found)
data_exp_filtered = data_exp.loc[data_exp['Unnamed: 0'].isin(DepMapId)].set_index('Unnamed: 0')
data_exp_filtered.index.name = 'DepMap_ID'
colnames = data_exp_filtered.columns
colnames_correct = [x.split()[0] for x in colnames]
data_exp_filtered.columns = colnames_correct

# Normalizing the expression (min-max) [ xscaled = (x - min(x)) / (max(x) - min(x)) ]
data_exp_norm = (data_exp_filtered-data_exp_filtered.min())/(data_exp_filtered.max()-data_exp_filtered.min())
data_exp_norm.index = \
    [data_info['stripped_cell_line_name'][data_info['DepMap_ID'] == x][0] for x in data_exp_norm.index]
data_exp_filtered.index = data_exp_norm.index

# Filter circ and neuro gene expression only (scaled)
clock_expression = data_exp_norm[circ_genes['Gene name']]
neuro_expression = data_exp_norm[neuro_genes['Gene name']]

# Clean CCLE drug sensitivity data
to_replace = ['-', '\\.', ' ', '/', ':', '\\(', '\\)']
for mark in to_replace:
    drug_sens['Primary Cell Line Name'] = drug_sens['Primary Cell Line Name'].replace(mark, '', regex=True)
drug_sens['Primary Cell Line Name'] = drug_sens['Primary Cell Line Name'].str.upper()
drug_sens['Compound'] = drug_sens['Compound'].replace('-', '', regex=True)

# Circadian score with cell line in upper cases
circ_score['Cellline'] = circ_score['Cellline'].str.upper()
circ_score.loc[2, 'Cellline'] = 'SKNBE2'

# Circadian feature values with cell line in upper cases
circ_values['cellline'] = circ_values['cellline'].str.upper()
circ_values.loc[2, 'cellline'] = 'SKNBE2'
circ_values = circ_values.set_index('cellline')
circ_values_norm = (circ_values-circ_values.min())/(circ_values.max()-circ_values.min())

# Reduce to drug columns of interest
drug_params = ['EC50 (uM)', 'IC50 (uM)', 'ActArea']  # not including Amax
drug_sens = drug_sens[['Primary Cell Line Name', 'Compound'] + drug_params]

# Find the available cell lines in CCLE drug data
drug_data = pd.DataFrame(columns=drug_sens.columns)
for line in celllines:
    if line in set(drug_sens['Primary Cell Line Name']):
        drug_data = pd.concat([drug_data, drug_sens.loc[drug_sens['Primary Cell Line Name'] == line]])

ccle_drugs_clean = drug_data['Compound'].unique()  # 23 drugs
remaining_celllines = drug_data['Primary Cell Line Name'].unique()  # 5 cell lines found

# We found 'SKNBE2', 'SKNSH', 'SKNAS', 'KELLY', 'CHP212'
# Others available: SKNFI, KPNSI9S, SIMA, SKNDZ, IMR32 (IMR5 and GIMEN missing, not as expected)
# In total only 10 'AUTONOMIC_GANGLIA' available

# For consistency throughout the project, we need at least 5 cell lines, no NaNs, and at least 3 unique values by drugs
drug_availability_dict = {measure: [] for measure in drug_params}
for measure in drug_params:
    drug_meas_here = drug_data[['Primary Cell Line Name', 'Compound', measure]].set_index(['Primary Cell Line Name'])
    for drug in ccle_drugs_clean:
        this_drug_here = drug_meas_here.loc[drug_meas_here['Compound'] == drug][measure]
        if len(this_drug_here) == 5 and ~this_drug_here.isna().any() and this_drug_here.nunique() >= 3:
            drug_availability_dict[measure].append(drug)

#####################
# ## 1) Clustering ##

# Clustering of circadian clock gene expression (columns) in neuroblastoma cell lines (rows).
# Clustermap plot (raw CCLE expression values, log2(TPM+1))
raw_clock_expression = data_exp_filtered[circ_genes['Gene name']]

min_val = np.min(np.ravel(raw_clock_expression))
max_val = np.max(np.ravel(raw_clock_expression))
center = np.median(np.ravel(raw_clock_expression))

fig = sns.clustermap(raw_clock_expression, cmap=plt.get_cmap('coolwarm'), figsize=(12, 6),
                     metric='euclidean', method='complete', dendrogram_ratio=(0.2, 0.2),
                     annot=False, vmin=min_val, center=center, vmax=max_val, linewidth=.003)
fig.ax_heatmap.set_ylabel(None)
# add frame lines around heatmap
fig.ax_heatmap.axhline(y=0, color='k', linewidth=2)
fig.ax_heatmap.axhline(y=raw_clock_expression.shape[0], color='k', linewidth=2)
fig.ax_heatmap.axvline(x=0, color='k', linewidth=2)
fig.ax_heatmap.axvline(x=raw_clock_expression.shape[1], color='k', linewidth=2)
fig.ax_col_dendrogram.set_title(f"Clustermap circadian gene expression\n(Raw values, log2(TPM+1))")
# Adapt color bar ticks
cbar = fig.ax_heatmap.collections[0].colorbar
cbar.ax.tick_params(labelsize=8)
cbar.ax.set_yticks([min_val, center, max_val], ['%.2f' % min_val, '%.2f' % center, '%.2f' % max_val])
# Save figure
for form in formats:
    plt.savefig(f'./{res_folder}/Figure_S4_a.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

######################
# ## 2) Correlation ##

# Spearman between core clock gene expression and drug sensitivity,including linear regression of significant hits

# ## Using the Bonferroni corrected p-value threshold
for measure in drug_params:
    # Actual computation
    drug_meas_here = drug_data[['Primary Cell Line Name', 'Compound', measure]].set_index(['Primary Cell Line Name'])
    corr_matrix_here = pd.DataFrame(index=drug_availability_dict[measure], columns=circ_genes['Gene name'])
    pval_matrix_here = pd.DataFrame(index=drug_availability_dict[measure], columns=circ_genes['Gene name'])
    for drug in corr_matrix_here.index:
        this_drug_here = drug_meas_here.loc[drug_meas_here['Compound'] == drug][measure]
        for gene in circ_genes['Gene name']:
            expression_here = clock_expression.loc[remaining_celllines, gene]
            corr_matrix_here.loc[drug, gene] = spearmanr(this_drug_here, expression_here).statistic
            pval_matrix_here.loc[drug, gene] = spearmanr(this_drug_here, expression_here).pvalue
    # Plot
    fig = sns.clustermap(corr_matrix_here.astype(float), cmap=plt.get_cmap('coolwarm'),
                         figsize=(corr_matrix_here.shape[1] / 1.5, corr_matrix_here.shape[0] / 1.5),
                         metric='euclidean', method='complete', dendrogram_ratio=(0.1, 0.1),
                         annot=False, vmin=corr_matrix_here.min().min(), center=0,
                         vmax=corr_matrix_here.max().max(), lw=.003)
    fig.ax_heatmap.set_ylabel(None)
    plt.setp(fig.ax_heatmap.get_yticklabels(), rotation=0)
    # Add frame lines around heatmap
    fig.ax_heatmap.axhline(y=0, color='k', linewidth=2)
    fig.ax_heatmap.axhline(y=corr_matrix_here.shape[0], color='k', linewidth=2)
    fig.ax_heatmap.axvline(x=0, color='k', linewidth=2)
    fig.ax_heatmap.axvline(x=corr_matrix_here.shape[1], color='k', linewidth=2)
    bonferroni_here = 0.05 / (len(drug_availability_dict[measure]) * len(circ_genes['Gene name']))
    fig.ax_col_dendrogram.set_title(
        f"Relationship between circadian gene expression and drug sensitivity\n"
        f"({measure.split()[0]}, Spearman r, {corr_matrix_here.shape[0]} drugs, "
        f"Bonferroni corrected p-value threshold = %.2e)" % bonferroni_here)
    correct_x = [x.get_text() for x in fig.ax_heatmap.get_xticklabels()]
    correct_y = [y.get_text() for y in fig.ax_heatmap.get_yticklabels()]
    pval_matrix_here_reordered = pval_matrix_here.reindex(index=correct_y, columns=correct_x)
    corr_matrix_here_reordered = corr_matrix_here.reindex(index=correct_y, columns=correct_x)
    for (i, j), z in np.ndenumerate(np.array(pval_matrix_here_reordered)):
        if z != 0 and z < 0.05:
            star_multiplier = int(('%.2e' % z)[-1]) - 1
            if star_multiplier > 4 or ('%.2e' % z)[-2:] == 10:
                star_multiplier = 4  # block if 4 or above
            t = '*' * star_multiplier
            fig.ax_heatmap.text(j + 0.5, i + 0.5, "%s\n%.1f" % (t, corr_matrix_here_reordered.iloc[i, j]),
                                ha="center", va="center", fontsize=10, weight="bold",
                                color='w' if np.abs(float(
                                    '%.1f' % corr_matrix_here_reordered.iloc[i, j])) >= 0.7 else 'k')
        if z != 0 and z < bonferroni_here:
            rect = patches.Rectangle((j, i), 1, 1, fill=False, edgecolor="black", lw=2)
            fig.ax_heatmap.add_patch(rect)
    # Adapt color bar ticks
    cbar = fig.ax_heatmap.collections[0].colorbar
    cbar.ax.tick_params(labelsize=8)
    cbar.ax.set_yticks([-1, 0, 1], [-1, 0, 1])
    # Save figure, correlations, and p values
    output_file = f'./{res_folder}/4_a_correlation_and_p_values.xlsx' if measure.split(' ')[0] == 'ActArea' else \
        f'./{res_folder}/S5_a_correlation_and_p_values.xlsx' if measure.split(' ')[0] == 'IC50' else \
        f'./{res_folder}/S5_b_correlation_and_p_values.xlsx'
    with pd.ExcelWriter(output_file) as writer:
        corr_matrix_here_reordered.to_excel(writer, sheet_name='correlation', index=False)
        pval_matrix_here_reordered.to_excel(writer, sheet_name='p-values', index=False)
    for form in formats:
        plt.savefig(f'./{res_folder}/Figure_4_a.{form}' if measure.split(' ')[0] == 'ActArea' else
                    f'./{res_folder}/Figure_S5_a.{form}' if measure.split(' ')[0] == 'IC50' else
                    f'./{res_folder}/Figure_S5_b.{form}', bbox_inches='tight', dpi=dpi)
    plt.close()


# Linear regression of the configuration that throw the most significant hits [(ActArea, Spearman)]
# - PER2-Topotecan, NR1D1-Paclitaxel, CLOCK-PF2341066, RORC-PLX4720 (spearman, 4 hits)

# - SPEARMAN hits
measure = 'ActArea'
combs_of_interest = {'PER2': 'Topotecan', 'NR1D1': 'Paclitaxel', 'CLOCK': 'PF2341066', 'RORC': 'PLX4720'}  # Gene : Drug

for gene, drug in combs_of_interest.items():
    drug_meas_here = drug_data[['Primary Cell Line Name', 'Compound', measure]].set_index(['Primary Cell Line Name'])
    this_drug_here = drug_meas_here.loc[drug_meas_here['Compound'] == drug][measure]
    this_drug_here.name = drug
    expression_here = clock_expression.loc[remaining_celllines, gene]
    df_reg = pd.concat([this_drug_here, expression_here], axis=1)
    df_hue = df_reg.reset_index()
    fig = sns.lmplot(x=gene, y=drug, data=df_hue, fit_reg=False, hue='index', scatter_kws={'alpha': 1})
    sns.regplot(x=gene, y=drug, data=df_reg, fit_reg=True, ci=95, scatter=False, truncate=True,
                line_kws={'color': 'k', 'lw': 1})
    sns.move_legend(fig, "upper right", bbox_to_anchor=(1.02, 0.58), frameon=True, title='NB cell lines')
    # Annotate stats
    stats = spearmanr(df_reg[drug], df_reg[gene])
    bonferroni_threshold = 0.05/(20 * 16)  # 20*16 we know this for act area
    plt.figtext(0.84, 0.62, f"R² = %.1f | r = %.1f\nBf$_{{adj ρ}}$ = %.1e" %
                (np.sqrt(np.abs(stats.statistic)), stats.statistic, stats.pvalue),
                ha="left", va="center", fontsize=10)
    plt.xlabel(f'{gene} gene expression')
    plt.ylabel(f'{measure} {drug}')
    plt.title('Linear regression')
    plt.tight_layout()
    # Save figure
    for form in formats:
        plt.savefig(f'./{res_folder}/Figure_4_c_{drug}.{form}', bbox_inches='tight', dpi=dpi)
    plt.close()

#######################
# ## 3) LDA analysis ##

# 3b, 3c, 3d, S4c input: neuro expression, target: clock score
# 3e, 3f, 3g S4e input: clock expression, target: clock score

combination = ['neuro_score', 'clock_score']
inputs = [neuro_expression, clock_expression]
targets = [circ_score, circ_score]

# STEP 1 ## Scatter plot and bar plot using all available features

# Scatter plot and bar plot
for comb, inp, target in zip(combination, inputs, targets):
    # Figure 3b and 3e
    expression_here = inp
    target_here = target.set_index('Cellline').reindex(inp.index)
    circ_score_here = target_here.loc[inp.index, 'mra_circadian']
    this_score_binarized = circ_score_here >= circ_score_here.median()
    subtypes_here = pd.DataFrame(list(circ_score_here.index), columns=['subtype'], index=circ_score_here.index)
    loocv_acc, _, _, _, _, _, not_correct = \
        LDA_loocv(data=expression_here, y=pd.concat([expression_here, this_score_binarized], axis=1),
                  target='mra_circadian')
    lda = LinearDiscriminantAnalysis(n_components=1)
    LDAs = lda.fit_transform(expression_here, this_score_binarized)
    LDAdf = pd.DataFrame(data=LDAs, columns=['LD1'], index=expression_here.index)
    _, _, ratio = get_LDA_metrics(this_score_binarized, LDAdf)
    # Now we can start the plotting
    plot_LDA(expression_here, pd.concat([expression_here, this_score_binarized, subtypes_here], axis=1),
             label='mra_circadian',
             title=f'LDA, median-binarized, circ. score, neuro expression\n'
                   f'ratio={"%.1f" % ratio}, LOOCV acc.={"%.2f" % loocv_acc}\n'
                   f'Not correct: {list(not_correct)}' if comb == 'neuro_score' else
                   f'LDA, median-binarized, circ. score, clock expression\n'
                   f'ratio={"%.1f" % ratio}, LOOCV acc.={"%.2f" % loocv_acc}\nNot correct: {list(not_correct)}',
             lda_output_target='mra_circadian', seed=seed, annot_subtype=True, annot_cellline=False,
             data_info=data_info, expand=(1.5, 2.2))
    # Save figure clock_drug
    for form in formats:
        plt.savefig(f'./{res_folder}/Figure_3_b.{form}' if comb == 'clock_score' else
                    f'./{res_folder}/Figure_3_e.{form}', bbox_inches='tight', dpi=dpi)
    plt.close()

# STEP 2 ## For 3c, 3d, 3f, 3g, we analyse the best gene expression combinations, plot their bcd-wcd ratios and best LDA
# including their LOOCV accuracy. 3f and 3g show top 4 genes, S4 bcde show top 1 - 6, S4b and S4d show loocv meta.
cells = ['Accuracy (%)', 'SKNSH', 'GIMEN', 'NGP', 'SKNAS', 'KELLY', 'SHSY5Y', 'SKNBE2', 'CHP212']
top_genes_circ = [1, 2, 3, 4, 5, 6, len(clock_expression.columns)]
top_genes_neuro = [1, 2, 3, 4, 5, 6, len(neuro_expression.columns)]

loocv_meta_circ = pd.DataFrame(index=cells, columns=top_genes_circ)
loocv_meta_neuro = pd.DataFrame(index=cells, columns=top_genes_neuro)
# Scatter plot and bar plot ranking on drug-by-drug basis
for comb, inp, target in zip(combination, inputs, targets):
    print(f'Processing LOOCV LDA in the {comb} combination...')
    expression_here = inp
    gene_combos = [comb for r in range(1, len(expression_here.columns) + 1) for comb in
                   itertools.combinations(expression_here.columns, r)]
    gene_combo_store = pd.DataFrame(index=['bcd-wcd_ratio'], columns=gene_combos)
    target_here = target.set_index('Cellline').reindex(inp.index)
    circ_score_here = target_here.loc[inp.index, 'mra_circadian']
    this_score_binarized = circ_score_here >= circ_score_here.median()
    subtypes_here = pd.DataFrame(list(circ_score_here.index), columns=['subtype'], index=circ_score_here.index)
    # Now we can start the plotting
    for g in gene_combos:
        expression_tmp = expression_here[list(g)]
        lda = LinearDiscriminantAnalysis(n_components=1)
        LDAs = lda.fit_transform(expression_tmp, this_score_binarized)
        LDAdf = pd.DataFrame(data=LDAs, columns=['LD1'], index=expression_tmp.index)
        _, _, ratio = get_LDA_metrics(this_score_binarized, LDAdf)
        gene_combo_store[g] = np.log10(ratio)
    # Plot all and highlight the best combination
    f, ax = plt.subplots()
    ax.scatter(np.arange(len(gene_combos)), gene_combo_store.values, s=20)
    ax.axhline(y=0, color='grey', linestyle='--', linewidth=1)  # marking the zero line
    ax.axhline(y=np.log10(2), color='grey', linestyle='-', linewidth=2)  # marking the ratio = 2 line
    combo_max_name = gene_combo_store.max(numeric_only=True).idxmax()
    ymax = gene_combo_store.max(numeric_only=True).max()
    # we are only interested in the best combinations up to the overall best one, len(combo_max_name)
    txts = []
    for gene_length in np.arange(1, len(combo_max_name) + 1):
        cols_here = [col for col in gene_combo_store.columns if len(col) == gene_length]
        combo_max_here = gene_combo_store[cols_here].max(numeric_only=True).idxmax()
        xmax_here = gene_combos.index(combo_max_here)
        ymax_here = gene_combo_store[cols_here].max(numeric_only=True).max()
        plt.scatter(xmax_here, ymax_here, marker='x', s=30, c='red')
        txts.append(plt.text(xmax_here, ymax_here, list(combo_max_here), ha='center', va='center'))
    adjust_text(txts, expand=(0, 2.5), force_text=(1.5, 1.5), only_move={'text': 'y-'},
                arrowprops=dict(arrowstyle='-', color='k'))
    plt.title(comb.replace('_', ' ~ ') +
              f' ({len(expression_here.columns)} genes, {len(gene_combos)} combinations)')
    plt.xlabel('# combination')
    plt.ylabel('BCD-WCD ratio [log 10]')
    for form in formats:
        plt.savefig(f'./{res_folder}/Figure_3_c.{form}' if comb == 'clock_score' else
                    f'./{res_folder}/Figure_3_f.{form}', bbox_inches='tight', dpi=dpi)
    plt.close()
    # Plot the best LDA with up to 6 genes, store also the 11/16 for the loocv meta plot
    gene_length_to_plot = 6
    for gene_length in np.append(np.arange(1, gene_length_to_plot + 1), len(expression_here.columns)):
        cols_here = [col for col in gene_combo_store.columns if len(col) == gene_length]
        combo_max_here = gene_combo_store[cols_here].max(numeric_only=True).idxmax()
        ymax_here = gene_combo_store[cols_here].max(numeric_only=True).max()
        loocv_acc, _, _, _, _, _, not_correct = \
            LDA_loocv(data=expression_here[list(combo_max_here)],
                      y=pd.concat([expression_here[list(combo_max_here)], this_score_binarized], axis=1),
                      target='mra_circadian')
        if comb == 'clock_score':
            loocv_meta_circ.loc['Accuracy (%)', gene_length] = loocv_acc * 100
            for cell in cells[1:]:
                loocv_meta_circ.loc[cell, gene_length] = 0 if cell in not_correct else 1
        if comb == 'neuro_score':
            loocv_meta_neuro.loc['Accuracy (%)', gene_length] = loocv_acc * 100
            for cell in cells[1:]:
                loocv_meta_neuro.loc[cell, gene_length] = 0 if cell in not_correct else 1
        if gene_length <= gene_length_to_plot:
            plot_LDA(expression_here[list(combo_max_here)],
                     pd.concat([expression_here[list(combo_max_here)], this_score_binarized, subtypes_here], axis=1),
                     label='mra_circadian',
                     title=f'LDA, median-binarized, circ. score, neuro expression'
                           f'\n{list(combo_max_here)}, ratio log10={"%.1f" % ymax_here}, '
                           f'LOOCV acc.={"%.2f" % loocv_acc}\n'
                           f'Not correct: {list(not_correct)}' if comb == 'neuro_score' else
                           f'LDA, median-binarized, circ. score, clock expression'
                           f'\n{list(combo_max_here)}, ratio log10={"%.1f" % ymax_here}, '
                           f'LOOCV acc.={"%.2f" % loocv_acc}\nNot correct: {list(not_correct)}',
                     lda_output_target='mra_circadian', seed=seed, annot_subtype=True, annot_cellline=False,
                     data_info=data_info, expand=(1.5, 2.2))
            for form in formats:
                plt.savefig(f'./{res_folder}/Figure_3_d.{form}' if comb == 'clock_score' and gene_length == 4 else
                            f'./{res_folder}/Figure_3_g.{form}' if comb == 'neuro_score' and gene_length == 4 else
                            f'./{res_folder}/Figure_S4_c_top_{gene_length}.{form}' if comb == 'clock_score' else
                            f'./{res_folder}/Figure_S4_e_top_{gene_length}.{form}', bbox_inches='tight', dpi=dpi)
            plt.close()
    print(f'Done with the {comb} combination.\n')

# plot the loocv meta heatmap
for meta, comb in zip([loocv_meta_circ, loocv_meta_neuro], ['clock_score', 'neuro_score']):
    mask = np.zeros_like(meta, dtype=bool)
    mask[0, :] = True
    cmap = sns.color_palette(["lightgrey", "cornflowerblue"])
    fig = sns.heatmap(meta.astype(float), mask=mask, cmap=cmap, cbar=False, linewidths=0.5, linecolor='black',
                      annot=False, fmt="", annot_kws={"fontsize": 12, "color": "black"})
    for i in range(meta.shape[1]):
        fig.text(i + 0.5, 0.5, f'{meta.iloc[0, i]:.1f}',
                 ha='center', va='center', color='black', fontsize=8)
    fig.xaxis.tick_top()
    fig.xaxis.set_label_position('top')
    plt.title('LOOCV - circadian genes' if comb == 'clock_score' else 'LOOCV - NB-specific genes')
    plt.tight_layout()
    for form in formats:
        plt.savefig(f'./{res_folder}/Figure_S4_b.{form}' if comb == 'clock_score' else
                    f'./{res_folder}/Figure_S4_d.{form}', bbox_inches='tight', dpi=dpi)
    plt.close()

##############################################################################################
# END OF SCRIPT FOR ANALYSIS PIPELINE FOR THE NEUROBLASTOMA PAPER COLLAB WITH CHARITE-BERLIN #
##############################################################################################
