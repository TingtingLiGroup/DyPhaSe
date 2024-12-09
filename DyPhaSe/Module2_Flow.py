import os
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import networkx as nx
import community as community_louvain
import random

parent_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def process_trend_clustering_with_r(input_file, output_dir, cluster_num=6, random_seed=123, 
                                    cutoff=0.85, r_home=None, r_user=None):
    """
    Perform trend clustering analysis on genes using the Mfuzz package in R after filtering genes
    based on their max DyPhaSe scores.

    Parameters:
        input_file (str): Path to the input DyPhaSe score result file (CSV format, e.g., DyPhaSe_scores.csv).
        output_dir (str): Directory to save the output results.
        cluster_num (int): Number of clusters for the analysis.
        random_seed (int): Random seed for reproducibility.
        cutoff (float): Threshold to filter out genes with DyPhaSe scores lower than this value. Default is 0.85.
        r_home (str, optional): Path to the R installation. If not provided, the system's default R installation is used.
                                Note that the required R packages (e.g., Mfuzz) must be installed in this environment.
        r_user (str, optional): Path to the R library containing the required packages.

    Returns:
        None: Results are saved directly in the specified output directory.
    """
    if r_home:
        os.environ["R_HOME"] = r_home
    if r_user:
        os.environ["R_USER"] = r_user

    pandas2ri.activate()

    r(f'''
    library(widgetTools, lib.loc="{r_user}")
    library(tkWidgets, lib.loc="{r_user}")
    ''')

    mfuzz = importr("Mfuzz", lib_loc=r_user)
    dplyr = importr("dplyr")
    biobase = importr("Biobase", lib_loc=r_user)

    r_code = f'''
    function(input_file, output_dir, cluster_num=6, random_seed=123, cutoff=0.85) {{
        library(Mfuzz)
        library(dplyr)

        df <- read.csv(input_file, row.names = 1, check.names = FALSE)
        df <- as.matrix(df)

        df <- df[apply(df, 1, function(row) any(row > {cutoff})), ]
        print(dim(df))

        mfuzz_class <- new('ExpressionSet', exprs = df)

        mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
        mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
        mfuzz_class <- filter.std(mfuzz_class, min.std = 0)

        mfuzz_class <- standardise(mfuzz_class)

        set.seed(random_seed)
        mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))

        cluster_centers <- mfuzz_cluster$centers
        cluster_centers_df <- as.data.frame(cluster_centers)
        rownames(cluster_centers_df) <- paste0("Cluster", 1:cluster_num)
        colnames(cluster_centers_df) <- colnames(df)
        write.csv(cluster_centers_df, file.path(output_dir, 'mfuzz-result-centers.csv'), row.names = TRUE)

        write.table(mfuzz_cluster$membership, file.path(output_dir, 'mfuzz-result-membership.csv'), sep = ',', col.names = NA, quote = FALSE)

        protein_cluster <- mfuzz_cluster$cluster
        protein_cluster <- cbind(df[names(protein_cluster), ], protein_cluster)
        write.table(protein_cluster, file.path(output_dir, 'mfuzz-result-cluster.csv'), sep = ',', col.names = NA, quote = FALSE)
    }}
    '''
    r_func = r(r_code)
    r_func(input_file, output_dir, cluster_num, random_seed, cutoff)




def plot_dyphase_clusters(centers, membership, target_ls=None, save_path=None, figsize=(5, 2)):
    """
    Plot trend curves for DyPhaSe clustering.

    Parameters:
        centers (pd.DataFrame): Cluster center data, where rows represent cluster IDs and columns represent time points.
        membership (pd.DataFrame): Membership data, where rows represent proteins and columns represent cluster IDs.
        target_ls (list, optional): List of target cluster IDs to be plotted. Defaults to all cluster IDs.
        save_path (str, optional): Path to save the image. If None, the plot will be displayed instead of saved.
        figsize (tuple, optional): Figure size, default is (5, 2).
    """
    if target_ls is None:
        target_ls = centers.index.tolist()
    
    centers = centers.loc[target_ls]
    membership = membership.loc[:, membership.columns.isin(target_ls)]
    time_points = centers.columns
    plt.figure(figsize=figsize, dpi=300)

    for cluster_id in centers.index:
        y = centers.loc[cluster_id].values
        max_membership_indices = membership.idxmax(axis=1) == cluster_id
        member_values = membership.loc[max_membership_indices, cluster_id].values
        x = np.arange(len(time_points))
        X_Y_Spline = make_interp_spline(x, y)
        X_ = np.linspace(x.min(), x.max(), 500)
        Y_ = X_Y_Spline(X_)
        std_deviation = np.std(member_values)
        lower_bound = y - std_deviation
        upper_bound = y + std_deviation
        lower_spline = make_interp_spline(x, lower_bound)(X_)
        upper_spline = make_interp_spline(x, upper_bound)(X_)

        plt.plot(X_, Y_, label=f'{cluster_id}', linewidth=2, alpha=0.8)
        plt.fill_between(X_, lower_spline, upper_spline, alpha=0.2)
        plt.scatter(x, y, color=plt.gca().lines[-1].get_color(), zorder=5, s=20)

    plt.xticks(ticks=np.arange(len(time_points)), labels=time_points, fontsize=12)
    plt.xlabel('', fontsize=12)
    plt.ylabel('Normalized\nDyPhaSe score', fontsize=12)
    plt.title('', fontsize=14)
    plt.yticks(ticks=[-1, 0, 1])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend(title="Clusters", loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    plt.grid(False)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    else:
        plt.show()



def filter_clusters_with_expression(
    gene_expression, 
    cell_types, 
    mfuzz_results, 
    cluster_conditions, 
    output_dir, 
    rmean_cutoff=0.2
):
    """
    Generate a table of average gene expression in each cell type and filter mfuzz clustering results 
    based on expression levels at DyPhaSe peak time. Removes genes with low expression at DyPhaSe peak time.

    Parameters:
        gene_expression (pd.DataFrame): Gene expression matrix with rows as cell names and columns as gene names.
        cell_types (pd.DataFrame): Cell type data with rows as cell names and a 'celltype' column indicating cell types.
        mfuzz_results (pd.DataFrame): mfuzz clustering results.
        cluster_conditions (dict): Filtering conditions for each cluster, e.g., {1: ['CT26', 'CT34'], 2: ['CT14', 'CT18']}.
        output_dir (str): Directory to save the results.
        rmean_cutoff (float, optional): Threshold for filtering based on rmean, default is 0.2.

    Returns:
        pd.DataFrame: Filtered clustering results.
    """
    os.makedirs(output_dir, exist_ok=True)

    gene_expression = gene_expression.T
    gene_expression = gene_expression.groupby(gene_expression.index).first()
    gene_expression = gene_expression[[col for col in gene_expression.columns]]

    cell_type_list = list(set(cell_types['celltype']))

    sta_data_expression = pd.DataFrame()

    for cell_type in cell_type_list:
        cell_indices = cell_types[cell_types['celltype'] == cell_type].index.tolist()
        cell_data = gene_expression[cell_indices]

        sta_data_expression[f'{cell_type}_mean'] = cell_data.apply(lambda x: np.mean(list(x)), axis=1)
        sta_data_expression[f'{cell_type}_rmean'] = (
            sta_data_expression[f'{cell_type}_mean']
            .rank(ascending=True, method='dense')
        )
        sta_data_expression[f'{cell_type}_rmean'] /= sta_data_expression[f'{cell_type}_rmean'].max()

    sta_data_expression.index.name = 'Gene'

    expression_path = os.path.join(output_dir, "sta_data_expression.csv")
    sta_data_expression.to_csv(expression_path, index=True)
    print(f"Intermediate result (sta_data_expression) saved to {expression_path}.")

    sta_data_expression = sta_data_expression[[col for col in sta_data_expression.columns if col.endswith('_rmean')]]
    filtered_df_mfuzz = pd.DataFrame()

    for cluster, conditions in cluster_conditions.items():
        sub_cluster = mfuzz_results[mfuzz_results['protein_cluster'] == cluster]

        condition = pd.Series([True] * len(sta_data_expression), index=sta_data_expression.index)
        for cond in conditions:
            condition &= sta_data_expression[f'{cond}_rmean'] > rmean_cutoff

        lowexp_ls = sta_data_expression[~condition].index.tolist()
        valid_indices = [idx for idx in sub_cluster.index if idx not in lowexp_ls]

        filtered_df_mfuzz = pd.concat([filtered_df_mfuzz, sub_cluster.loc[valid_indices]])

    final_path = os.path.join(output_dir, "filtered_clusters.csv")
    filtered_df_mfuzz.to_csv(final_path, index=True)
    print(f"Final filtered result saved to {final_path}.")

    return filtered_df_mfuzz



def analyze_ppi_community(df_mfuzz, species, clusters, node_num_cutoff=10, output_path=None, seed=123):
    """
    Analyze protein-protein interaction (PPI) networks and perform community detection for specified clusters. 
    Integrates the results of all clusters into a single file.

    Parameters:
        df_mfuzz (pd.DataFrame): DataFrame containing clustering information.
        species (str): Species information ("Human" or "Mouse").
        clusters (list): List of cluster IDs to analyze.
        node_num_cutoff (int, optional): Minimum node count threshold for filtering communities. Default is 10.
        output_path (str, optional): Path to save the results. If None, results are not saved.
        seed (int, optional): Random seed for reproducibility. Default is 123.

    Returns:
        pd.DataFrame: Combined results of all analyzed clusters.
    """
    np.random.seed(seed)
    random.seed(seed)

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

    all_results = []

    for cluster_id in clusters:
        protein_list = df_mfuzz[df_mfuzz['protein_cluster'] == cluster_id].index.tolist()

        ppi_nw = pd.read_csv(species_files[species]["ppi"])
        ppi_nw.columns = ['A', 'B']
        treat_proteins = set(protein_list)
        ppi_nw = ppi_nw[ppi_nw['A'].isin(treat_proteins) & ppi_nw['B'].isin(treat_proteins)]
        ppi_nw[['A', 'B']] = np.sort(ppi_nw[['A', 'B']], axis=1)
        ppi_nw = ppi_nw.drop_duplicates(subset=['A', 'B'])

        G = nx.from_pandas_edgelist(ppi_nw, 'A', 'B')
        partition = community_louvain.best_partition(G)
        community_sizes = {comm_id: sum(1 for node in partition if partition[node] == comm_id) for comm_id in set(partition.values())}
        final_partition = {node: comm_id for node, comm_id in partition.items() if community_sizes[comm_id] >= node_num_cutoff}
        filtered_G = G.subgraph(final_partition.keys()).copy()

        unique_community_ids = sorted(set(final_partition.values()))
        new_id_mapping = {old_id: new_id for new_id, old_id in enumerate(unique_community_ids, 1)}
        renamed_partition = {node: new_id_mapping[comm_id] for node, comm_id in final_partition.items()}

        for node, community_id in renamed_partition.items():
            filtered_G.nodes[node]['community'] = community_id

        community_clustering_dict = {}
        for community_id, nodes in pd.DataFrame.from_dict(renamed_partition, orient='index').groupby(0):
            community_nodes = nodes.index.tolist()
            if len(community_nodes) > 1:
                subgraph = filtered_G.subgraph(community_nodes)
                community_clustering_dict[community_id] = nx.average_clustering(subgraph)

        valid_communities = set(community_clustering_dict.keys())
        community_dict = {
            new_id: [node for node in renamed_partition if renamed_partition[node] == new_id]
            for new_id in valid_communities
        }

        protein_community_pairs = [
            (protein, comm_id)
            for comm_id, proteins in community_dict.items()
            for protein in proteins
        ]
        protein_community_df = pd.DataFrame(protein_community_pairs, columns=['Protein', 'Community'])

        degree_centrality = nx.degree_centrality(filtered_G)
        betweenness_centrality = nx.betweenness_centrality(filtered_G)
        closeness_centrality = nx.closeness_centrality(filtered_G)
        eigenvector_centrality = nx.eigenvector_centrality(filtered_G)

        protein_community_df['Degree Centrality'] = protein_community_df['Protein'].map(degree_centrality)
        protein_community_df['Betweenness Centrality'] = protein_community_df['Protein'].map(betweenness_centrality)
        protein_community_df['Closeness Centrality'] = protein_community_df['Protein'].map(closeness_centrality)
        protein_community_df['Eigenvector Centrality'] = protein_community_df['Protein'].map(eigenvector_centrality)
        protein_community_df['SaPS'] = protein_community_df['Protein'].apply(lambda x: gp_dic.get(x, np.nan)).map(saps_dic)
        protein_community_df['Community Clustering'] = protein_community_df['Community'].map(community_clustering_dict)
        protein_community_df['Cluster'] = cluster_id

        all_results.append(protein_community_df)

    all_results_df = pd.concat(all_results, ignore_index=True)
    all_results_df['name'] = 'Cluster' + all_results_df['Cluster'].astype(str) + '-c' + all_results_df['Community'].astype(str)

    if output_path:
        all_results_df.to_csv(f'{output_path}/target_cluster_ppi_community_node_info.csv', index=False)
        print(f"All results saved to {output_path}")

    return all_results_df




def generate_cytoscape_files(name_ls, protein_community_df, output_dir):
    """
    Generate subnet files for Cytoscape visualization.

    Parameters:
        name_ls (list): List of community names.
        protein_community_df (pd.DataFrame): DataFrame containing protein and community names.
        output_dir (str): Directory to save the output files.

    Returns:
        None
    """

    for name in name_ls:
        protein_list = protein_community_df[protein_community_df['name'] == name]['Protein'].tolist()
        human_ppi_nw = pd.read_csv(os.path.join(parent_folder, "data/10090.string.ppi.genename.csv"))
        human_ppi_nw.columns = ['A', 'B']
        treat_proteins = set(protein_list)
        human_ppi_nw = human_ppi_nw[human_ppi_nw['A'].isin(treat_proteins) & human_ppi_nw['B'].isin(treat_proteins)]
        human_ppi_nw[['A', 'B']] = np.sort(human_ppi_nw[['A', 'B']], axis=1)
        human_ppi_nw = human_ppi_nw.drop_duplicates(subset=['A', 'B'])
        output_file_path = f"{output_dir}/{name}-ppi-pairs.csv"
        human_ppi_nw.to_csv(output_file_path, index=False)
        print(f"Saved: {output_file_path}")
