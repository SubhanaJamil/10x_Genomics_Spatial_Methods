

<h1>Analysis of Visium Fluorescence Data Using Squidpy</h1>

<div class="section">
    <h2>1. Introduction</h2>
    <p>
        Spatial transcriptomics allows the study of gene expression while preserving tissue architecture.
        In this analysis, we used Visium fluorescence data from a mouse brain section to integrate
        image-based features with gene expression clusters.
    </p>
    <p>The workflow combines:</p>
    <ul>
        <li>Image processing and segmentation</li>
        <li>Feature extraction (summary, histogram, texture)</li>
        <li>Clustering based on image-derived features</li>
    </ul>
    <p>
        The aim is to determine whether morphological features from images provide complementary
        biological insights beyond gene expression.
    </p>
</div>

<div class="section">
    <h2>2. Importing Libraries and Loading Data</h2>
    <p>
        This section initializes all required libraries for spatial and image analysis. The dataset is loaded
        using Squidpy’s built-in Visium fluorescence example. The AnnData object contains gene expression
        data, while the ImageContainer stores the tissue image. Logging ensures correct setup.
    </p>
    <pre><code>%matplotlib inline
import pandas as pd
import matplotlib.pyplot as plt

!pip install anndata
!pip install scanpy
!pip install squidpy

import anndata as ad
import scanpy as sc
import squidpy as sq

sc.logging.print_header()
print(f"squidpy=={sq.__version__}")

img = sq.datasets.visium_fluo_image_crop()
adata = sq.datasets.visium_fluo_adata_crop()</code></pre>
</div>

<div class="section">
    <h2>3. Spatial Visualization of Clusters</h2>
    <p>
        Spatial gene expression clusters are visualized using spatial_scatter. Each spot corresponds
        to a tissue location and is colored by cluster identity.
    </p>
    <pre><code>sq.pl.spatial_scatter(adata, color="cluster")</code></pre>
</div>


<img width="642" height="491" alt="Spatial Visualization of Clusters" src="https://github.com/user-attachments/assets/ea3beaec-9c50-4a6e-8155-377324493883" />


<div class="section">
    <h2>4. Visualization of Fluorescence Channels</h2>
    <p>
        The dataset includes DAPI, NEUN, and GFAP channels. These help identify DNA, neurons,
        and glial cells based on fluorescence intensity.
    </p>
    <pre><code>img.show(channelwise=True)</code></pre>
</div>

<img width="790" height="289" alt="Visualization of Fluorescence Channels" src="https://github.com/user-attachments/assets/71ea47cd-cdfd-42ce-9ae2-da1f9087bbd8" />


<div class="section">
    <h2>5. Image Preprocessing and Segmentation</h2>
    <p>
        The image is smoothed to reduce noise, followed by watershed segmentation to detect nuclei.
        Each segmented region represents a cell or nucleus.
    </p>
    <pre><code>sq.im.process(
    img=img,
    layer="image",
    method="smooth",
)

sq.im.segment(img=img, layer="image_smooth", method="watershed", channel=0, chunks=1000)

fig, ax = plt.subplots(1, 2)
img_crop = img.crop_corner(2000, 2000, size=500)
img_crop.show(layer="image", channel=0, ax=ax[0])
img_crop.show(layer="segmented_watershed", channel=0, ax=ax[1])</code></pre>
</div>

<div class="section">
    <h2>6. Extraction of Segmentation Features</h2>
    <p>
        Segmentation features quantify morphology such as cell count and intensity, helping estimate
        cell density and identify cell types.
    </p>
    <pre><code>features_kwargs = {"segmentation": {"label_layer": "segmented_watershed"}}

sq.im.calculate_image_features(
    adata,
    img,
    features="segmentation",
    layer="image",
    key_added="features_segmentation",
    n_jobs=1,
    features_kwargs=features_kwargs,
)

sq.pl.spatial_scatter(
    sq.pl.extract(adata, "features_segmentation"),
    color=[
        "segmentation_label",
        "cluster",
        "segmentation_ch-0_mean_intensity_mean",
        "segmentation_ch-1_mean_intensity_mean",
    ],
    frameon=False,
    ncols=2)
</code></pre>
</div>
<p align='center'> 
<img width="515" height="267" alt="Segmentation Feature Extraction" src="https://github.com/user-attachments/assets/a35d7629-1a9c-4cd7-b333-8e68a200b059" />
</p>

<div class="section">
    <h2>7. Extraction of Multi-scale Image Features</h2>
    <p>
        Multi-scale features (summary, histogram, texture) capture both local and global patterns
        in tissue morphology.
    </p>
    <pre><code>params = {
    "features_orig": {
        "features": ["summary", "texture", "histogram"],
        "scale": 1.0,
        "mask_circle": True,
    },
    "features_context": {"features": ["summary", "histogram"], "scale": 1.0},
    "features_lowres": {"features": ["summary", "histogram"], "scale": 0.25},
}

for feature_name, cur_params in params.items():
    sq.im.calculate_image_features(
        adata, img, layer="image", key_added=feature_name, n_jobs=1, **cur_params
    )

adata.obsm["features"] = pd.concat(
    [adata.obsm[f] for f in params.keys()], axis="columns"
)

adata.obsm["features"].columns = ad.utils.make_index_unique(
    adata.obsm["features"].columns
)</code></pre>
</div>

<img width="1088" height="829" alt="Multi_scale Feature Extraction" src="https://github.com/user-attachments/assets/ca24f968-b942-493c-ac95-1b53354f4074" />

<div class="section">
    <h2>8. Clustering Based on Image Features</h2>
    <p>
        PCA and Leiden clustering are applied to image features to generate clusters based on morphology.
    </p>
    <pre><code>def cluster_features(features: pd.DataFrame, like=None):
    if like is not None:
        features = features.filter(like=like)

   adata = ad.AnnData(features)
    sc.pp.scale(adata)
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)

   return adata.obs["leiden"]

adata.obs["features_summary_cluster"] = cluster_features(
    adata.obsm["features"], like="summary"
)
adata.obs["features_histogram_cluster"] = cluster_features(
    adata.obsm["features"], like="histogram"
)
adata.obs["features_texture_cluster"] = cluster_features(
    adata.obsm["features"], like="texture"
)</code></pre>
</div>

<div class="section">
    <h2>9. Visualization and Comparison of Clusters</h2>
    <h3>Explanation</h3>
    <p>
        Image-based clusters are compared with gene expression clusters to assess biological relevance.
    </p>

  <h3>Code</h3>
    <pre><code>sc.set_figure_params(facecolor="white", figsize=(8, 8))

sq.pl.spatial_scatter(
    adata,
    color=[
        "features_summary_cluster",
        "features_histogram_cluster",
        "features_texture_cluster",
        "cluster",
    ],
    ncols=3,
)
</code></pre>
</div>

<img width="3435" height="2205" alt="Cluster Comparison" src="https://github.com/user-attachments/assets/c26657be-4b4b-44d8-9116-9c6f057f86a2" />


<div class="section">
    <h2>10. Results</h2>
    <ul>
        <li>Segmentation identified nuclei and enabled estimation of cell density.</li>
        <li>Hippocampus regions showed higher cell density.</li>
        <li>NEUN → neuron-rich regions</li>
        <li>GFAP → glial-rich regions</li>
        <li>Image-based clustering showed higher spatial resolution.</li>
        <li>Captured fine anatomical structures like cortical layers.</li>
    </ul>
</div>

<div class="section">
    <h2>11. Conclusion</h2>
    <p>
        Integrating image features with gene expression enhances spatial transcriptomics analysis.
        Gene-based clustering provides broad classification, while image features provide fine-grained insights.
    </p>
    <ul>
        <li>Image features complement gene expression</li>
        <li>Segmentation improves interpretation</li>
        <li>Multi-scale features capture tissue complexity</li>
    </ul>
    <p>
        Squidpy-based image analysis is a powerful tool for extracting biological meaning from spatial datasets.
    </p>
</div>

</body>
</html>
