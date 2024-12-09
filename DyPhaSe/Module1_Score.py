import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


parent_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def calculate_dyphase_score(gene_expression, cell_types, species, output_dir):
    """
    Calculate the stability score, differential stability score, and DyPhaSe score for each gene
    across different cell types, and save the results.

    Parameters:
        gene_expression (pd.DataFrame): Gene expression matrix, where rows are cell names and columns are gene names.
        cell_types (pd.DataFrame): Cell type information for each cell, indexed by cell name with a 'celltype' column specifying the cell type.
        species (str): Species information for the single-cell data's gene names, either 'Human' or 'Mouse'.
        output_dir (str): Directory where the results will be saved.

    Returns:
        tuple: A tuple containing three DataFrames:
            - df_s: Stability scores for each gene.
            - df_ds: Differential stability scores for each gene.
            - df_rpsds: DyPhaSe scores for each gene.
    """
    os.makedirs(output_dir, exist_ok=True)

    species_files = {
        "Human": {
            "gene_protein": os.path.join(parent_folder, "data/human_uniprot_gene_ensembl.csv"),
            "phase_scores": os.path.join(parent_folder, "data/human_LLPS_score.csv"),
            "ppi": os.path.join(parent_folder, "data/9606.string.ppi.genename.csv")
        },
        "Mouse": {
            "gene_protein": os.path.join(parent_folder, "data/mouse_uniprot_gene_ensembl.csv"),
            "phase_scores": os.path.join(parent_folder, "data/mouse_LLPS_score.csv"),
            "ppi": os.path.join(parent_folder, "data/10090.string.ppi.genename.csv")
        }
    }
    
    if species not in species_files:
        raise ValueError("Invalid species name. Choose 'Human' or 'Mouse'.")

    gene_protein_data = pd.read_csv(species_files[species]["gene_protein"])
    gp_dic = dict(zip(gene_protein_data['Gene'], gene_protein_data['Uniprot_protein']))

    phase_scores_data = pd.read_csv(species_files[species]["phase_scores"])
    saps_dic = dict(zip(phase_scores_data['Protein'], phase_scores_data['rnk_SaPS']))

    if not gene_expression.index.equals(cell_types.index):
        raise ValueError("Indices mismatch: Rows of gene_expression and cell_types must match exactly.")

    valid_genes = [
        gene for gene in gene_expression.columns 
        if gene in gp_dic and gp_dic[gene] in saps_dic
    ]
    gene_expression = gene_expression[valid_genes]
    gene_expression = gene_expression.loc[:, (gene_expression != 0).any(axis=0)]

    unique_cell_types = cell_types['celltype'].unique()
    df_s = pd.DataFrame()

    for cell_type in unique_cell_types:
        cells_of_type = cell_types[cell_types['celltype'] == cell_type].index
        subset_expression = gene_expression.loc[cells_of_type]
        subset_expression = subset_expression.loc[:, (subset_expression != 0).any(axis=0)]

        zero_proportion = (subset_expression == 0).sum(axis=0) / subset_expression.shape[0]
        mean_expression = subset_expression.mean(axis=0)
        sorted_genes = mean_expression.sort_values().index
        std_dev = subset_expression.loc[:, sorted_genes].std(axis=0)
        rollm_std = std_dev.rolling(window=50, center=False, min_periods=25).median()
        stddm = std_dev - rollm_std

        rank_zero_proportion = zero_proportion.rank()
        rank_stddm = stddm.rank()
        stability_score = 1 - (rank_zero_proportion + rank_stddm) / (2 * subset_expression.shape[1])

        df_s[cell_type] = stability_score
        df_s.dropna(inplace=True)

    df_ds = df_s.apply(lambda col: col - df_s.drop(columns=col.name).mean(axis=1), axis=0)
    df_ds['saps'] = [saps_dic[gp_dic[gene]] for gene in df_ds.index]
    df_psds = df_ds.drop(columns='saps').multiply(df_ds['saps'], axis=0)
    df_rpsds = df_psds.rank(axis=0) / df_psds.shape[0]

    df_s.to_csv(os.path.join(output_dir, "stability_scores.csv"))
    df_ds.to_csv(os.path.join(output_dir, "differential_stability_scores.csv"))
    df_rpsds.to_csv(os.path.join(output_dir, "DyPhaSe_scores.csv"))

    return df_s, df_ds, df_rpsds



def plot_dyphase_change(sta_data, gene_list, figsize=(6, 2), ymin=None, ymax=None, save_path=None, yticks=None):
    """
    Plot line graphs showing the change of DyPhaSe scores for specific genes across cell types. 
    The plot can be saved to a specified path or displayed directly.

    Parameters:
        - sta_data: pd.DataFrame, the DyPhaSe_scores.csv data.
        - gene_list: list, a list of genes to be plotted.
        - figsize (tuple, optional): Figure size, default is (6, 2).
        - ymin (float, optional): Minimum value for the y-axis. Default is None.
        - ymax (float, optional): Maximum value for the y-axis. Default is None.
        - save_path (str, optional): File path to save the plot. If None, the plot is displayed.
        - yticks (list, optional): Custom ticks for the y-axis. Default is None.
    """
    plt.figure(figsize=figsize, dpi=400)
        
    color_list = ['blue', 'orange', 'green', 'red', 'purple',
                  'brown', 'pink', 'grey', 'lightblue', 'lightcoral',
                  'lightgreen', 'lightpink', 'lightsalmon', 'lightseagreen',
                  'lightskyblue', 'lightslategray', 'lightsteelblue', 'lightyellow', 'lightcyan']
        
    sta_data.columns = [x.split('_')[0] for x in sta_data.columns.tolist()]

    j = 0
    for i, row in sta_data.loc[gene_list].iterrows():
        plt.plot(row.index, row.values, label=i, marker='o', 
                 markerfacecolor=color_list[j], linewidth=2.5, linestyle='dashed')
        j += 1

    plt.title('')
    plt.xlabel('')
    plt.ylabel('DyPhaSe Score', fontsize=12, fontproperties='Arial')
    plt.xlim(-0.5, len(sta_data.columns.tolist()) - 0.5)
    
    if ymin is None:
        ymin = 0
    if ymax is None:
        ymax = 1
    plt.ylim(ymin, ymax)
    
    if yticks is not None:
        plt.yticks(yticks, fontsize=12, fontproperties='Arial')
    else:
        plt.yticks(fontsize=12, fontproperties='Arial')

    legend_font = {'family': 'Arial', 'size': 12}
    legend = plt.legend(title='', bbox_to_anchor=(1.01, 1), prop=legend_font)
    
    plt.xticks(fontsize=12, fontproperties='Arial')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=400)
    else:
        plt.show()
