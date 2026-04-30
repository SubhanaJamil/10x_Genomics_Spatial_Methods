<h1>Spatial Transcriptomics Analysis using 10x Genomics Visium and Scanpy</h1>
<h2>1. Introduction</h2>
    <p>
        Spatial transcriptomics is a cutting-edge technique that enables the study of gene expression 
        within the spatial context of intact tissues. Unlike traditional RNA sequencing, which loses 
        positional information, spatial approaches preserve tissue architecture and reveal how gene 
        expression varies across different regions.
    </p>

   <p>
        The Visium platform developed by 10x Genomics allows high-throughput spatial profiling by 
        capturing mRNA from tissue sections placed on barcoded arrays. For computational analysis, 
        we used Scanpy, which provides efficient tools for preprocessing, clustering, and visualization 
        of transcriptomics data.
    </p>

  <h2>2. Dataset Description</h2>

  <p><strong>Dataset:</strong> Human Lymph Node (Visium)</p>
    <p><strong>Source:</strong> 10x Genomics</p>
 <ul>
        <li><strong>Initial Cells (Spots):</strong> 4035</li>
        <li><strong>Genes:</strong> 36,601</li>
    </ul>
 <p>
        The dataset contains gene expression values mapped to spatial coordinates along with 
        histological tissue images, enabling integrated spatial and molecular analysis.
 </p>
 <h2>3. Methodology</h2>

<h3>Step 1: Installation and Imports</h3>
<p>
This step ensures that all required libraries such as Scanpy, Pandas, Matplotlib, and Seaborn 
are installed and imported. These libraries are essential for handling biological data, 
performing statistical analysis, and generating visualizations.
</p>

<pre><code>!pip install scanpy

from __future__ import annotations
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
</code></pre>
<p>All dependencies were successfully installed and imported without errors.</p>

<h3>Step 2: Scanpy Settings</h3>
<p>
Scanpy settings were configured to control verbosity and figure aesthetics. This ensures 
reproducibility and consistent visualization output across the analysis pipeline.
</p>

<pre><code>sc.settings.verbosity = 3
sc.set_figure_params(facecolor="white", figsize=(8, 8))
</code></pre>

<p>Environment configured for detailed logging and standardized plots.</p>


<h3>Step 3: Data Loading</h3>
<p>
The Visium spatial dataset was loaded using Scanpy’s built-in dataset function. The gene 
names were made unique to avoid duplication issues.
</p>

<pre><code>adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata.var_names_make_unique()
</code></pre>

<p><strong>Result:</strong></p>
<ul>
    <li>Dataset loaded successfully</li>
    <li>AnnData object created with 4035 cells × 36,601 genes</li>
</ul>



<h3>Step 4: Quality Control Metrics</h3>
<p>
Quality control metrics were calculated to assess the quality of cells. Mitochondrial genes 
were identified, and their proportion was computed, which is an important indicator of cell health.
</p>

<pre><code>adata.var["mt"] = adata.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    inplace=True
)
</code></pre>

<p> QC metrics such as total counts, gene counts, and mitochondrial percentage were successfully added to the dataset.</p>


<h3>Step 5: QC Visualization</h3>
<p>
Histograms were generated to visualize the distribution of total counts and gene counts. 
These plots help in identifying thresholds for filtering low-quality cells.
</p>

<pre><code>fig, axs = plt.subplots(1, 4, figsize=(15, 4))

sns.histplot(adata.obs["total_counts"], ax=axs[0])
sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] &lt; 10000], bins=40, ax=axs[1])
sns.histplot(adata.obs["n_genes_by_counts"], bins=60, ax=axs[2])
sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] &lt; 4000], bins=60, ax=axs[3])

plt.show()
</code></pre>
<img width="2038" height="609" alt="QC Histogram" src="https://github.com/user-attachments/assets/a17298fb-77d1-4afc-8792-f0a6ec8bfdf9" />

<p align='center'> Distributions revealed appropriate thresholds for filtering based on counts and gene expression.</p>

<h3>Step 6: Filtering</h3>
<p>
Cells and genes that do not meet quality criteria were removed. This step eliminates noise 
and improves downstream analysis accuracy.
</p>

<pre><code>sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)

adata = adata[adata.obs["pct_counts_mt"] &lt; 20].copy()

sc.pp.filter_genes(adata, min_cells=10)
</code></pre>

<p><strong>Results:</strong></p>
<p>During quality control, low-quality cells and uninformative genes were removed to improve the reliability of downstream analysis. 44 cells with fewer than 5000 counts were filtered out, as these likely represent low-quality or damaged spots with insufficient RNA content. 
130 cells with more than 35000 counts were removed, as very high counts may indicate doublets (multiple cells captured in one spot), which can distort biological signals. 
  After mitochondrial filtering, the dataset retained 3861 high-quality cells, ensuring that stressed or dying cells were excluded.
  16916 genes expressed in fewer than 10 cells were also removed, since such lowly expressed genes provide little statistical power and mostly add noise.</p>
<table align='center' border="1" cellpadding="5">
    <tr>
        <th>Filtering Step</th>
        <th>Cells Remaining</th>
    </tr>
    <tr>
        <td>Initial</td>
        <td>4035</td>
    </tr>
    <tr>
        <td>After filtering</td>
        <td>3861</td>
    </tr>
</table>

<p>Overall, this filtering step reduced technical noise, removed unreliable cells and genes, and improved the quality of the dataset for accurate clustering and spatial analysis.</p>

<h3>Step 7: Normalization</h3>
<p>
Normalization was performed to remove technical variability between cells. Log transformation 
was applied to stabilize variance, and highly variable genes were identified.
</p>

<pre><code>sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat",
    n_top_genes=2000
)
</code></pre>

<p><strong>Result:</strong></p>
<ul>
    <li>Data normalized successfully</li>
    <li>2000 highly variable genes selected</li>
</ul>


<h3>Step 8: Dimensionality Reduction</h3>
<p>
Dimensionality reduction techniques such as PCA and UMAP were applied to simplify the 
high-dimensional data while preserving important structure.
</p>

<pre><code>sc.pp.pca(adata)
sc.pp.neighbors(adata)

sc.tl.umap(adata)
</code></pre>

<p>Low-dimensional representation of data generated for visualization and clustering.</p>



<h3>Step 9: Clustering</h3>
<p>
The Leiden clustering algorithm was used to group cells based on gene expression similarity. 
This helps identify distinct biological populations.
</p>

<pre><code>sc.tl.leiden(
    adata,
    key_added="clusters",
    flavor="igraph",
    directed=False,
    n_iterations=2
)
</code></pre>

<p>Total 10 Clusters were found.</p>

<h3>Step 10: UMAP Visualization</h3>
<p>
UMAP plots were generated to visualize clusters and QC metrics in reduced dimensional space.
</p>

<pre><code>sc.pl.umap(
    adata,
    color=["total_counts", "n_genes_by_counts", "clusters"],
    wspace=0.4
)
</code></pre>
<img width="4332" height="1591" alt="UMAP visualization" src="https://github.com/user-attachments/assets/36ddfbb7-61c8-4a9f-879b-4b8699a4ee04" />

<p align='center'> Distinct clusters were observed, indicating successful separation of cell populations.</p>


<h3>Step 11: Spatial Visualization</h3>
<p>
Spatial plots were generated to map gene expression and clusters onto the tissue image. 
This step integrates molecular and spatial information.
</p>

<pre><code>sc.pl.spatial(
    adata,
    img_key="hires",
    color=["total_counts", "n_genes_by_counts"]
)
  
</code></pre>
<img width="2269" height="1091" alt="spatial Visualization" src="https://github.com/user-attachments/assets/aff53b72-7478-4376-8946-62c51e316c21" />
<p align='center'> Clusters showed clear spatial organization within the tissue, confirming biological relevance.</p>
<p>
The function <code>scanpy.pl.spatial()</code> supports several important parameters that allow 
fine control over spatial visualization of gene expression data.
</p>

<p><strong>Key Parameters:</strong></p>
<ul>
    <li><code>img_key</code>: Key where the image is stored in <code>adata.uns</code></li>
    <li><code>crop_coord</code>: Coordinates for cropping the image (left, right, top, bottom)</li>
    <li><code>alpha_img</code>: Controls transparency of the tissue image</li>
    <li><code>bw</code>: Converts the image into grayscale (black & white)</li>
</ul>

<p>
Additionally, the <code>size</code> parameter behaves differently in spatial plots. 
Instead of defining absolute size, it acts as a scaling factor for spot sizes.
</p>

<hr>

<p>
Before spatial mapping, clustering was performed in gene expression space and visualized using UMAP. 
By projecting these clusters back onto spatial coordinates, we can better understand tissue organization 
and potential inter-cellular communication.
</p>

<pre><code>sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    size=1.5
)
</code></pre>
<p align='center'>
<img width="586" height="449" alt="S_visualization2" src="https://github.com/user-attachments/assets/e7ed2122-7089-454a-bed8-dfa55ff51140" />
</p>

<p align='center'>Spatial visualization revealed that clustered cells are not randomly distributed, 
but instead show structured organization across the tissue section, supporting biologically meaningful spatial patterns.</p>


<h3>Step 12: Marker Gene Identification</h3>

<p>
Marker genes were identified for each cluster using statistical testing. This helps determine 
which genes are highly expressed in specific cell populations.
</p>

<pre><code>sc.tl.rank_genes_groups(
    adata,
    "clusters",
    method="t-test"
)

sc.pl.rank_genes_groups_heatmap(
    adata,
    groups="9",
    n_genes=10,
    groupby="clusters"
)
</code></pre>


<p> Marker genes were successfully identified for all clusters, highlighting distinct gene expression signatures across cell populations.</p>



<h3>Step 13: Marker Gene Spatial Visualization</h3>

<p>
Identified marker genes were visualized in spatial coordinates to understand their distribution 
within the tissue and to explore biological relevance.
</p>

<pre><code>sc.pl.spatial(
    adata,
    img_key="hires",
    color=["clusters", "CR2"]
)
  </code></pre>
  <img width="2294" height="1091" alt="Marker gene visualization" src="https://github.com/user-attachments/assets/d9c3c84a-1a25-4c91-8a60-43c65fb102c5" />
<hr>
<pre><code>
sc.pl.spatial(
    adata,
    img_key="hires",
    color=["COL1A2", "SYPL1"],
    alpha=0.7
)
</code></pre>
<img width="2239" height="1091" alt="Marker gene visualization 2" src="https://github.com/user-attachments/assets/a222fc9d-6851-44bc-9afa-a44f4008bc40" />


<p> Genes such as CR2, COL1A2, and SYPL1 showed distinct spatial expression patterns, confirming cluster-specific biological activity.</p>


<h2>4. Conclusion</h2>
<p>
This study successfully implemented a complete spatial transcriptomics pipeline using Scanpy 
and Visium data from 10x Genomics.
</p>

<p><strong>Key Findings:</strong></p>
<ul>
    <li>Effective quality control and filtering</li>
    <li>Identification of 10 distinct clusters</li>
    <li>Spatially meaningful gene expression patterns</li>
</ul>

<p>
The integration of spatial and transcriptomic data provides valuable insights into tissue 
organization and cellular heterogeneity, demonstrating the power of spatial omics in modern 
bioinformatics.
</p>
