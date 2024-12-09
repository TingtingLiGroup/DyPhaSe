# Introduction to DyPhaSe
DyPhaSe (Dynamic Phase Separation Predictor) is a modular computational tool designed to model and predict phase separation dynamics based on gene expression variability. It provides a flexible framework for identifying phase-separating proteins and characterizing condensate behavior across diverse cellular states. By integrating protein-protein interaction (PPI) networks, DyPhaSe enables precise prediction of dynamic condensate flows within specific spatiotemporal contexts, shedding light on their functional roles in various biological processes.

This repository contains the source code for DyPhaSe, empowering researchers to predict dynamic phase separation events in their specific biological conditions of interest.

---
# Install
```bash
# Step 1: Create a new conda virtual environment with Python 3.9
conda create --name DyPhaSe python=3.9

# Step 2: Activate the new environment
conda activate DyPhaSe

# Step 3: Install the required libraries
pip install pandas==2.1.1 numpy==1.26.1 seaborn==0.11.0 matplotlib==3.8.1 rpy2==3.5.17 scipy==1.11.3 networkx==3.2.1 louvain==0.8.1

# Step 4: Clone the DyPhaSe repository from GitHub
git clone https://github.com/TingtingLiGroup/DyPhaSe.git

# Step 5: Navigate into the cloned repository directory
cd DyPhaSe

# Step 6: Install the package using setup.py
python setup.py install
```

---
# Function Documentation
## Module1: Predict Dynamic PS Ability of Single Protein
![alt text](./data/DyphaSe-figure%20(11).png)
In Module1, we provide the calculation function `calculate_dyphase_score` and the plotting function `plot_dyphase_change`.

### `calculate_dyphase_score`
#### Description
  The `calculate_dyphase_score` function calculates the predicted DyPhaSe score for each gene in different cell types based on the user-submitted single-cell gene expression matrix and cell type annotation files. A higher DyPhaSe score corresponds to an increased likelihood of phase separation under specific cellular conditions.

#### Parameters
+ **gene_expression** (*pd.DataFrame*): Gene expression matrix, where rows are cell names and columns are gene names. Make share the matrix has been normalized, such as log(TPM+1), log(CPM+1), etc.
+ **cell_types** (*pd.DataFrame*): Cell type information for each cell, indexed by cell name with a 'celltype' column specifying the cell type.
+ **species** (str): Species information for the single-cell data's gene names, either **'Human'** or **'Mouse'**.
+ **output_dir** (str): Directory where the results will be saved.

#### Returns
+ **df_s** (*pd.DataFrame*): Stability scores for each gene.
+ **df_ds** (*pd.DataFrame*): Differential stability scores for each gene.
+ **df_rpsds** (*pd.DataFrame*): DyPhaSe scores for each gene.

#### Example
```python
gene_expression = pd.read_csv('../data/gene_expression.csv', index_col=0)
cell_types = pd.read_csv('../data/cell_types.csv', index_col=0)
calculate_dyphase_score(gene_expression, 
                        cell_types, 
                        species='Mouse', output_dir='./DyPhaSe-package/out/')
```
---
### `plot_dyphase_change`
#### Description
  The `plot_dyphase_change` function plots line graphs showing the change of DyPhaSe scores for a given gene list across cell types.
  The plot can be saved to a specified path or displayed directly.

#### Parameters
+ **sta_data** (*pd.DataFrame*): The DyPhaSe_scores.csv data.
+ **gene_list** (*list*): A list of genes to be plotted.
+ **figsize** (*tuple, optional*): Figure size. Default is (6, 2).
+ **ymin** (*float, optional*): Minimum value for the y-axis. Default is None.
+ **ymax** (*float, optional*): Maximum value for the y-axis. Default is None.
+ **save_path** (*str, optional*): File path to save the plot. If None, the plot is displayed.
+ **yticks** (*list, optional*): Custom ticks for the y-axis. Default is None.

#### Example
```python
sta_data = pd.read_csv(f'../out/DyPhaSe_scores.csv',index_col=0)
gene_list = [ 'Pcnt', 'Dynll1','Dcdc2a', 'Tubgcp2', 'Tubb4a', 'Mapk8ip3','Ckap5']
plot_dyphase_change(sta_data, gene_list)
```

---
## Module2: Predict Condensate Flow Dynamically
![alt text](<./data/DyphaSe-figure (10).png>)
In Module2, we provide the calculation function `process_trend_clustering_with_r`, `filter_clusters_with_expression`, `analyze_ppi_community`, the plotting function `plot_dyphase_clusters` and the optional function `generate_cytoscape_files`.

### `process_trend_clustering_with_r`
#### Description
  The `process_trend_clustering_with_r` function performs trend clustering analysis on genes using the Mfuzz package in R after filtering genes based on their max DyPhaSe scores.
  Users need to make sure they have an R conda environment with R >= 4.4.0 and install required R packages: Mfuzz, dplyr, Biobase in order to run the function directly in python via rpy2.

#### Parameters
+ **input_file** (*str*): Path to the input DyPhaSe score result file (CSV format, e.g., DyPhaSe_scores.csv).
+ **output_dir** (*str*): Directory to save the output results.
+ **cluster_num** (*int*): Number of clusters for the analysis.
+ **random_seed** (*int*): Random seed for reproducibility.
+ **cutoff** (*float*): Threshold to filter out genes with DyPhaSe scores lower than this value. Default is 0.85.
+ **r_home** (*str, optional*): Path to the R installation. If not provided, the system's default R installation is used. Note that the required R packages (e.g., Mfuzz) must be installed in this environment.
+ **r_user** (*str, optional*): Path to the R library containing the required packages.

#### Returns
+ **'mfuzz-result-centers.csv'** (*pd.DataFrame*): Contains the cluster centers from the Mfuzz clustering, representing the central gene expression values for each cluster.
+ **'mfuzz-result-membership.csv'** (*pd.DataFrame*): Records the cluster assignment for each gene, indicating which cluster each gene belongs to.
+ **'mfuzz-result-cluster.csv'** (*pd.DataFrame*): Combines the gene expression data with cluster assignments, showing both the expression values and the cluster number for each gene.

#### Example
```python
process_trend_with_r(
    input_file="../out/DyPhaSe_scores.csv",
    output_dir="../out/",
    cluster_num=12,
    random_seed=123,
    cutoff=0.85, 
    r_home="/home/wyq/.conda/envs/R442/lib/R",  
    r_user="/home/wyq/.conda/envs/R442/lib/R/library" 
)
```
---
### `plot_dyphase_clusters`
#### Description
  The `plot_dyphase_clusters` function plots curves for each DyPhaSe clustering pattern.
  The plot can be saved to a specified path or displayed directly.

#### Parameters
- **centers** (*pd.DataFrame*): Cluster center data, where rows represent cluster IDs and columns represent time points.
- **membership** (*pd.DataFrame*): Membership data, where rows represent proteins and columns represent cluster IDs.
- **target_ls** (*list, optional*): List of target cluster IDs to be plotted. Defaults to all cluster IDs.
- **save_path** (*str, optional*): Path to save the image. If None, the plot will be displayed instead of saved.
- **figsize** (*tuple, optional*): Figure size, default is (5, 2).

#### Example
```python
centers = pd.read_csv('../out/mfuzz-result-centers.csv', index_col=0)
membership = pd.read_csv('../out/mfuzz-result-membership.csv', index_col=0)
plot_dyphase_clusters(centers, membership, target_ls=['Cluster1','Cluster9','Cluster10','Cluster11','Cluster12'])
```
---
### `filter_clusters_with_expression`
#### Description
  The `filter_clusters_with_expression` function generates a table of average gene expression in each cell type and filter mfuzz clustering results based on expression levels at DyPhaSe peak time. Removes genes with low expression at DyPhaSe peak time.

#### Parameters
- **gene_expression** (*pd.DataFrame*): Gene expression matrix with rows as cell names and columns as gene names.
- **cell_types** (*pd.DataFrame*): Cell type data with rows as cell names and a 'celltype' column indicating cell types.
- **mfuzz_results** (*pd.DataFrame*): Mfuzz clustering results.
- **cluster_conditions** (*dict*): Filtering conditions for each cluster, e.g., {1: ['CT26', 'CT34'], 2: ['CT14', 'CT18']}.
- **output_dir** (*str*): Directory to save the results.
- **rmean_cutoff** (*float, optional*): Threshold for filtering based on rmean, default is 0.2.

#### Returns
- **filtered_df_mfuzz** (*pd.DataFrame*): Filtered clustering results, containing only genes that meet the specified expression criteria at DyPhaSe peak time.

#### Example
```python
gene_expression = pd.read_csv('../data/gene_expression.csv', index_col=0)
cell_types = pd.read_csv('../data/cell_types.csv', index_col=0)
mfuzz_results = pd.read_csv('../out/mfuzz-result-cluster.csv',index_col=0)

cluster_conditions = {
    1: ['CT26', 'CT30'],
    2: ['CT26', 'CT34'],
    3: ['CT14', 'CT18', 'CT26'],
    4: ['CT18', 'CT26', 'CT30'],
    5: ['CT22', 'CT34'],
    6: ['CT14', 'CT22', 'CT30'],
    7: ['CT14', 'CT22'],
    8: ['CT14', 'CT18', 'CT30'],
    9: ['CT30', 'CT34'],
    10: ['CT22', 'CT26'],
    11: ['CT18', 'CT22'],
    12: ['CT14', 'CT18', 'CT34']
}

filtered_results = filter_clusters_with_expression(
    gene_expression=gene_expression,
    cell_types=cell_types,
    mfuzz_results=mfuzz_results,
    cluster_conditions=cluster_conditions,
    output_dir="../out/",
    rmean_cutoff=0.2
)
```

---
### `analyze_ppi_community`
#### Description
  The `analyze_ppi_community` function analyzes protein-protein interaction (PPI) networks and perform community detection by Louvain algorithm for genes in specified clusters.

#### Parameters
- **df_mfuzz** (*pd.DataFrame*): DataFrame containing clustering information.
- **species** (*str*): Species information ("Human" or "Mouse").
- **clusters** (*list*): List of cluster IDs to analyze.
- **node_num_cutoff** (*int, optional*): Minimum node count threshold for filtering communities. Default is 10.
- **output_path** (*str, optional*): Path to save the results. If None, results are not saved.
- **seed** (*int, optional*): Random seed for reproducibility. Default is 123.

#### Returns
- **all_results_df** (*pd.DataFrame*): Combined results of all analyzed clusters, including community detection and network centrality metrics for each protein in the specified clusters.

#### Example
```python
df_mfuzz = pd.read_csv('../out/filtered_clusters.csv',index_col=0)
combined_results = analyze_ppi_community(
    df_mfuzz=df_mfuzz,
    species="Mouse",
    clusters=[1,9,10,11,12],
    node_num_cutoff=10,
    output_path="../out/",
    seed=123
)
```
---
### `generate_cytoscape_files`
#### Description
  The `generate_cytoscape_files` generates a file containing interaction pairs of selected community for Cytoscape visualization.

#### Parameters
- **name_ls** (*list*): List of community names.
- **protein_community_df** (*pd.DataFrame*): DataFrame containing protein and community names.
- **output_dir** (*str*): Directory to save the output files.

#### Returns
- **'ppi-pairs.csv'** (*pd.DataFrame*): The CSV files containing protein-protein interaction pairs for each specified community in the output directory.

#### Example
```python
name_ls = ['Cluster11-c1', 'Cluster10-c1']
protein_community_df = pd.read_csv("../out/target_cluster_ppi_community_node_info.csv")
output_dir = "../out/"
generate_cytoscape_files(name_ls, protein_community_df, output_dir)
```