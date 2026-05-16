<div class="section">
<h1>Visium H&E Spatial Transcriptomics Analysis (Squidpy)</h1>
<p class="note">
This project demonstrates how to analyze spatial transcriptomics data using <b>Squidpy</b> and <b>Scanpy</b>.
We work with a mouse brain Visium H&E dataset and perform image feature extraction, clustering, spatial statistics, and ligand-receptor analysis.
</p>
</div>

<!-- 1 -->
<div class="section">
<h2>1. Install Required Libraries</h2>
<p class="note">
We install required bioinformatics libraries including Scanpy for single-cell analysis and Squidpy for spatial transcriptomics.
</p>

<pre><code>%matplotlib inline

!pip install anndata scanpy==1.11.1 squidpy leidenalg hdbscan==0.8.29 pandas==2.2.2
</code></pre>

<p class="note">
  This ensures all dependencies for spatial analysis and clustering are available.
</p>
</div>

<!-- 2 -->
<div class="section">
<h2>2. Import Libraries</h2>
<p class="note">
We import core Python libraries for data handling and biological analysis workflows.
</p>

<pre><code>import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc
import squidpy as sq
</code></pre>

<p class="note">
NumPy and Pandas handle data, while Scanpy and Squidpy handle single-cell and spatial data.
</p>
</div>

<!-- 3 -->
<div class="section">
<h2>3. Load Visium Dataset</h2>
<p class="note">
We load a preprocessed Visium H&E mouse brain dataset along with its histology image.
</p>

<pre><code>img = sq.datasets.visium_hne_image()
adata = sq.datasets.visium_hne_adata()
</code></pre>

<p class="note">
<code>adata</code> contains gene expression per spot, and <code>img</code> contains the tissue image.
</p>
</div>

<!-- 4 -->
<div class="section">
<h2>4. Spatial Visualization</h2>
<p class="note">
We visualize spatial clusters directly on the tissue section to understand anatomical structure.
</p>

<pre><code>sq.pl.spatial_scatter(
    adata,
    color="cluster"
)
</code></pre>

<p class="note">
This shows how different brain regions cluster spatially based on gene expression.
</p>
</div>
<img width="651" height="221" alt="spatial cluster" src="https://github.com/user-attachments/assets/bea0d17c-5c6f-4207-95d8-5d9989098aaf" />


<!-- 5 -->
<div class="section">
<h2>5. Image Feature Extraction</h2>
<p class="note">
We extract morphological features from histology images at multiple scales.
</p>

<pre><code>for scale in [1.0, 2.0]:

    feature_name = f"features_summary_scale{scale}"

    sq.im.calculate_image_features(
        adata,
        img.compute(),
        features="summary",
        key_added=feature_name,
        n_jobs=4,
        scale=scale,
    )
</code></pre>

<p class="note">
Larger scale captures broader tissue context, while smaller scale captures local structure.
</p>
</div>

<!-- 6 -->
<div class="section">
<h2>6. Combine Image Features</h2>
<p class="note">
We merge extracted features from multiple scales into one dataset for downstream analysis.
</p>

<pre><code>adata.obsm["features"] = pd.concat(
    [
        adata.obsm[f]
        for f in adata.obsm.keys()
        if "features_summary" in f
    ],
    axis="columns",
)

adata.obsm["features"].columns = ad.utils.make_index_unique(
    adata.obsm["features"].columns
)
</code></pre>

<p class="note">
This creates a unified feature matrix combining all image-derived information.
</p>
</div>

<!-- 7 -->
<div class="section">
<h2>7. Clustering Image Features</h2>
<p class="note">
We cluster image-derived features using PCA + Leiden clustering to find morphological patterns.
</p>

<pre><code>def cluster_features(features: pd.DataFrame, like=None):

    if like is not None:
        features = features.filter(like=like)

    adata_tmp = ad.AnnData(features)

    sc.pp.scale(adata_tmp)

    sc.pp.pca(adata_tmp, n_comps=min(10, features.shape[1] - 1))

    sc.pp.neighbors(adata_tmp)

    sc.tl.leiden(adata_tmp)

    return adata_tmp.obs["leiden"]
</code></pre>

<p class="note">
PCA reduces dimensionality, and Leiden groups similar tissue structures.
</p>
</div>
<img width="1469" height="431" alt="2nd" src="https://github.com/user-attachments/assets/1b886ab1-0333-40ce-befc-769e81a3ba9d" />

<!-- 8 -->
<div class="section">
<h2>8. Assign Feature Clusters</h2>
<p class="note">
We apply clustering results to each spatial spot.
</p>

<pre><code>adata.obs["features_cluster"] = cluster_features(
    adata.obsm["features"],
    like="summary"
)
</code></pre>

<p class="note">
Each spot now has both gene-based and image-based cluster labels.
</p>
</div>

<!-- 9 -->
<div class="section">
<h2>9. Compare Gene vs Image Clusters</h2>
<p class="note">
We compare clustering from gene expression vs tissue morphology.
</p>

<pre><code>sq.pl.spatial_scatter(
    adata,
    color=["features_cluster", "cluster"]
)
</code></pre>

<p class="note">
This helps identify agreement or mismatch between gene and image patterns.
</p>
</div>

<!-- 10 -->
<div class="section">
<h2>10. Spatial Graph Construction</h2>
<p class="note">
We build a spatial neighbor graph where nearby tissue spots are connected.
</p>

<pre><code>sq.gr.spatial_neighbors(adata)
</code></pre>

<p class="note">
This graph is required for spatial statistics and neighborhood analysis.
</p>
</div>

<!-- 11 -->
<div class="section">
<h2>11. Neighborhood Enrichment</h2>
<p class="note">
We test whether certain clusters tend to be neighbors more than expected by chance.
</p>

<pre><code>sq.gr.nhood_enrichment(
    adata,
    cluster_key="cluster"
)
</code></pre>

<p class="note">
 Helps identify biologically meaningful spatial interactions.
</p>
</div>

<!-- 12 -->
<div class="section">
<h2>12. Co-occurrence Analysis</h2>
<p class="note">
We compute how often clusters appear near each other in space.
</p>

<pre><code>sq.gr.co_occurrence(
    adata,
    cluster_key="cluster"
)
</code></pre>

<pre><code>sq.pl.co_occurrence(
    adata,
    cluster_key="cluster",
    clusters="Hippocampus",
    figsize=(8, 4),
)
</code></pre>

<p class="note">
  Shows spatial dependency between brain regions.
</p>
</div>
<img width="811" height="411" alt="3rd cooccurenec" src="https://github.com/user-attachments/assets/ee5eae93-bb16-46e8-bed8-3c1a6c60b38c" />

<!-- 13 -->
<div class="section">
<h2>13. Ligand-Receptor Interaction</h2>
<p class="note">
We predict cell-cell communication signals using ligand-receptor pairs.
</p>

<pre><code>sq.gr.ligrec(
    adata,
    n_perms=100,
    cluster_key="cluster",
)
</code></pre>

<pre><code>sq.pl.ligrec(
    adata,
    cluster_key="cluster",
    source_groups="Hippocampus",
    target_groups=[
        "Pyramidal_layer",
        "Pyramidal_layer_dentate_gyrus",
    ],
    means_range=(3, np.inf),
    alpha=1e-4,
    swap_axes=True,
)
</code></pre>

<p class="note">
 Identifies possible signaling between brain regions.
</p>
</div>

<!-- 14 -->
<div class="section">
<h2>14. Spatial Gene Patterns (Moran’s I)</h2>
<p class="note">
We identify genes that show spatially structured expression patterns.
</p>

<pre><code>genes = adata[:, adata.var.highly_variable].var_names.values[:1000]

sq.gr.spatial_autocorr(
    adata,
    mode="moran",
    genes=genes,
    n_perms=100,
    n_jobs=1,
)
</code></pre>

<p class="note">
Moran’s I detects genes that are clustered in space rather than randomly distributed.
</p>
</div>

<!-- 15 -->
<div class="section">
<h2>15. Visualize Spatial Gene Expression</h2>
<p class="note">
We visualize expression of spatially important genes on tissue.
</p>

<pre><code>sq.pl.spatial_scatter(
    adata,
    color=["Olfm1", "Plp1", "Itpka", "cluster"]
)
</code></pre>

<p class="note">
  Shows how specific genes map to brain regions.
</p>
</div>
<img width="2586" height="431" alt="spatial expression of selected 4th" src="https://github.com/user-attachments/assets/169f3bab-c030-4a38-9667-d48a9c7ad951" />


</body>
</html>
