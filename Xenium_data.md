

<h1>Analysis of Xenium Spatial Transcriptomics Data Using Squidpy and SpatialData</h1>


<h2>1. Introduction</h2>

<p>Spatial transcriptomics allows the study of gene expression while preserving tissue architecture. In this analysis, a <b>Xenium Human Kidney dataset</b> was used to explore spatial gene expression patterns using <b>SpatialData</b>, <b>Scanpy</b>, and <b>Squidpy</b>.</p>

<p>The workflow includes:</p>

<ul>
<li>Dataset download and loading</li>
<li>Quality control analysis</li>
<li>Filtering and normalization</li>
<li>Dimensionality reduction and clustering</li>
<li>Spatial statistics (co-occurrence, centrality, enrichment)</li>
<li>Gene expression visualization</li>
<li>Interactive spatial exploration</li>
</ul>

<p>The goal is to understand <b>cellular organization and spatial gene regulation in tissue microenvironment</b>.</p>


<h2>2. Environment Setup and Library Installation</h2>


<p>This section installs and imports all required libraries for spatial transcriptomics analysis. It ensures compatibility between SpatialData, Squidpy, and Scanpy. Additional tools like Matplotlib and Seaborn are used for visualization. Napari is installed for interactive spatial exploration.</p>


<pre>
get_ipython().system('pip install spatialdata spatialdata-io squidpy')
import spatialdata as sd
from spatialdata_io import xenium
import spatialdata_plot # Import spatialdata_plot to enable .pl accessor

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq
</pre>


<h2>3. Dataset Download and Loading</h2>

<p>The Xenium dataset is downloaded from the 10x Genomics repository, extracted, and loaded using the SpatialData Xenium reader. The object is stored as <b>sdata</b>, which contains images, shapes, and transcript data.</p>


<pre>
!mkdir -p Xenium
!wget -O Xenium/data.zip "https://cf.10xgenomics.com/samples/xenium/4.0.0/Xenium_V1_Protein_Human_Kidney_tiny/Xenium_V1_Protein_Human_Kidney_tiny_outs.zip"
!unzip Xenium/data.zip -d Xenium/
from spatialdata_io import xenium

sdata = xenium("./Xenium", cells_as_circles=True)
print(sdata)
</pre>



<h2>4. Zarr Conversion and Re-loading</h2>

<p>The dataset is converted into Zarr format for efficient storage and reloaded for analysis. This ensures faster access and scalability for large spatial datasets.</p>


<pre>
import zarr
zarr.config.set({'array.rectilinear_chunks': True})
sdata.write("./Xenium.zarr", overwrite=True)

from spatialdata_io import xenium

sdata = xenium("./Xenium", cells_as_circles=True)
sdata.write("Xenium.zarr", overwrite=True)
</pre>


<h2>5. Extracting AnnData Object</h2>

<p>The AnnData object contains gene expression counts, cell metadata, and spatial coordinates. It is extracted from the SpatialData object for downstream analysis.</p>


<pre>
adata = sdata.tables["table"]
adata
adata.obs
</pre>

<h2>Cell-Level Spatial Transcriptomics Summary</h2>

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; width: 100%;">
    <thead style="background-color: #f2f2f2;">
        <tr>
            <th>Cell ID</th>
            <th>Transcript Counts</th>
            <th>Total Counts</th>
            <th>Cell Area</th>
            <th>Nucleus Area</th>
            <th>Nucleus Count</th>
            <th>Segmentation Method</th>
            <th>Region</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td>aadmbfof-1</td>
            <td>29</td>
            <td>29</td>
            <td>43.89</td>
            <td>18.42</td>
            <td>1</td>
            <td>ATP1A1+CD45+E-Cad boundary</td>
            <td>cell_circles</td>
        </tr>
        <tr>
            <td>aageapbo-1</td>
            <td>46</td>
            <td>46</td>
            <td>42.22</td>
            <td>23.12</td>
            <td>1</td>
            <td>ATP1A1+CD45+E-Cad boundary</td>
            <td>cell_circles</td>
        </tr>
        <tr>
            <td>aakefffb-1</td>
            <td>27</td>
            <td>27</td>
            <td>18.24</td>
            <td>12.73</td>
            <td>1</td>
            <td>ATP1A1+CD45+E-Cad boundary</td>
            <td>cell_circles</td>
        </tr>
        <tr>
            <td>abkjennb-1</td>
            <td>29</td>
            <td>29</td>
            <td>41.04</td>
            <td>28.49</td>
            <td>1</td>
            <td>ATP1A1+CD45+E-Cad boundary</td>
            <td>cell_circles</td>
        </tr>
        <tr>
            <td>acdcmlfl-1</td>
            <td>40</td>
            <td>40</td>
            <td>22.08</td>
            <td>13.63</td>
            <td>1</td>
            <td>ATP1A1+CD45+E-Cad boundary</td>
            <td>cell_circles</td>
        </tr>
    </tbody>
</table>


<h2>6. QC Preparation for Spatial Analysis</h2>

<p>Spatial metadata is initialized to ensure compatibility with Squidpy plotting functions. Dummy spatial keys are added to prevent visualization errors during spatial plotting.</p>


<pre>
# Initialize adata.uns['spatial'] for squidpy plotting
if "spatial" not in adata.uns:
    adata.uns["spatial"] = {}

if "spatial" not in adata.uns["spatial"]:
    adata.uns["spatial"]["spatial"] = {}

adata.uns["spatial"]["spatial"]["images"] = {}
adata.uns["spatial"]["spatial"]["scalefactors"] = {
    "hires": 1.0,
    "lowres": 1.0,
    "spot_diameter_fullres": 1.0
}
adata.uns["spatial"]["spatial"]["metadata"] = {
    "obsm_key": "spatial"
}
</pre>

<h2>Spatial Coordinates of Cells</h2>

<table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; width: 100%;">
    <thead style="background-color: #f2f2f2;">
        <tr>
            <th>Cell Index</th>
            <th>X Coordinate</th>
            <th>Y Coordinate</th>
        </tr>
    </thead>
    <tbody>
        <tr><td>1</td><td>64.64</td><td>613.75</td></tr>
        <tr><td>2</td><td>107.50</td><td>621.21</td></tr>
        <tr><td>3</td><td>123.23</td><td>630.40</td></tr>
        <tr><td>4</td><td>116.32</td><td>627.98</td></tr>
        <tr><td>5</td><td>122.88</td><td>626.17</td></tr>
        <tr><td>6</td><td>119.04</td><td>622.41</td></tr>
        <tr><td>7</td><td>105.66</td><td>617.04</td></tr>
        <tr><td>8</td><td>97.07</td><td>626.52</td></tr>
        <tr><td>9</td><td>103.02</td><td>628.43</td></tr>
        <tr><td>10</td><td>100.68</td><td>614.05</td></tr>
        <tr><td>...</td><td>...</td><td>...</td></tr>
        <tr><td>358</td><td>65.03</td><td>36.72</td></tr>
    </tbody>
</table>


<h2>7. Quality Control Metrics</h2>

<p>QC metrics are calculated to assess dataset quality. The percentage of control probes and decoding errors indicates technical noise levels.</p>


<pre>
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)

cprobes = (
    adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
cwords = (
    adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")
</pre>

<h2>Quality Control Results</h2>

<p><b>Result:</b> During preprocessing of the Xenium spatial transcriptomics dataset, quality control was performed to assess background noise and decoding accuracy. Two key metrics were evaluated:</p>

<ul>
    <li><b>Negative DNA probe count (%):</b> 0.0116%</li>
    <li><b>Negative decoding count (%):</b> 0.0058%</li>
</ul>

<p>
These extremely low percentages indicate minimal background signal and very low decoding errors. This suggests that the dataset has high technical quality, with reliable probe hybridization and accurate barcode decoding.
</p>

<p>
Overall, the QC results confirm that the data is suitable for downstream spatial gene expression analysis.
</p>
<h2>8. QC Visualization</h2>

<p>Histograms show distributions of transcript counts, gene counts, cell area, and nucleus ratio. These help identify outliers and define filtering thresholds.</p>



<pre>
fig, axs = plt.subplots(1, 4, figsize=(15, 4))

axs[0].set_title("Total transcripts per cell")
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])

axs[1].set_title("Unique transcripts per cell")
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1])

axs[2].set_title("Area of segmented cells")
sns.histplot(adata.obs["cell_area"], kde=False, ax=axs[2])

axs[3].set_title("Nucleus ratio")
sns.histplot(adata.obs["nucleus_area"] / adata.obs["cell_area"], kde=False, ax=axs[3])
</pre>


<h2>9. Filtering and Preprocessing</h2>

<p>Low-quality cells and genes are filtered. Data is then normalized, log-transformed, and prepared for dimensionality reduction.</p>


<pre>
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_genes(adata, min_cells=5)

!pip install leidenalg
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
</pre>


<h2>10. Dimensionality Reduction and Clustering</h2>

<p>PCA reduces dimensionality, while UMAP visualizes data structure. Leiden clustering identifies cell populations based on transcriptomic similarity.</p>


<pre>
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
</pre>



<h2>11. UMAP Visualization</h2>



<pre>
sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
)
</pre>



<h2>12. Spatial Clustering Visualization</h2>



<pre>
sq.pl.spatial_scatter(
    adata,
    library_id="spatial",
    shape=None,
    color=[
        "leiden",
    ],
    wspace=0.4,
    img=False,
    size=2
)
</pre>


<h2>13. Spatial Graph and Centrality Analysis</h2>



<pre>
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

sq.gr.centrality_scores(adata, cluster_key="leiden")

sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(16, 5))
</pre>



<h2>14. Subsampling and Co-occurrence Analysis</h2>


<pre>
sdata.tables["subsample"] = sc.pp.subsample(adata, fraction=0.5, copy=True)

adata_subsample = sdata.tables["subsample"]

if 'spatial' not in adata_subsample.uns or adata_subsample.uns['spatial'] is None:
    adata_subsample.uns['spatial'] = adata.uns['spatial'].copy()

sq.gr.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
)

sq.pl.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
    clusters=["1", "2"],
    figsize=(10, 10),
)

sq.pl.spatial_scatter(
    adata_subsample,
    color="leiden",
    shape=None,
    size=2,
    library_id="spatial"
)
</pre>

<h2>15. Neighborhood Enrichment Analysis</h2>


<pre>
sq.gr.nhood_enrichment(adata, cluster_key="leiden")

fig, ax = plt.subplots(1, 2, figsize=(13, 7))
sq.pl.nhood_enrichment(
    adata,
    cluster_key="leiden",
    figsize=(8, 8),
    title="Neighborhood enrichment adata",
    ax=ax[0],
)
sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2, ax=ax[1])
</pre>

<h2>16. Spatial Autocorrelation (Moran’s I)</h2>


<pre>
sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)

sq.gr.spatial_autocorr(
    adata_subsample,
    mode="moran",
    n_perms=100,
    n_jobs=1,
)

adata_subsample.uns["moranI"].head(10)
</pre>

<h2>17. Gene Expression Visualization</h2>


<pre>
genes = list(adata.var_names[:2])
print(genes)

sq.pl.spatial_scatter(
    adata_subsample,
    library_id="spatial",
    color=genes,
    shape=None,
    size=2,
    img=False,
)
</pre>

<h2>18. Morphology-based Visualization</h2>

<pre>
sdata.pl.render_images("morphology_focus").pl.render_shapes(
    "cell_circles",
    table_name="table",
    use_raw=False,
    color="ACKR1"
).pl.show()

sdata.pl.render_images("morphology_focus").pl.render_shapes(
    "cell_circles",
    table_name="table",
    use_raw=False,
    color="ACTA2"
).pl.show()
</pre>


<h2>19. Interactive Visualization</h2>


<pre>
!pip install napari-spatialdata PyQt5
!Xvfb :99 -ac & display=:99

from napari_spatialdata import Interactive

Interactive(sdata)
</pre>


<h2>20. Results</h2>

<ul>
<li>Xenium dataset successfully loaded and processed</li>
<li>QC metrics confirmed high-quality transcript data</li>
<li>Leiden clustering identified distinct cell populations</li>
<li>Centrality scores highlighted important spatial hubs</li>
<li>Co-occurrence revealed cluster interactions</li>
<li>Moran’s I detected spatially variable genes</li>
<li>Morphology plots showed clear gene localization patterns</li>
<li>Interactive visualization enabled deeper exploration</li>
</ul>


<h2>21. Conclusion</h2>

<p>This study demonstrates an end-to-end spatial transcriptomics workflow using Squidpy, Scanpy, and SpatialData.</p>

<p>Key findings:</p>

<ul>
<li>Strong spatial organization of gene expression</li>
<li>Clear clustering of cell populations</li>
<li>Significant spatial gene patterns identified</li>
<li>Integration of morphology and transcriptomics improves biological insight</li>
</ul>

<p>Overall, Xenium spatial analysis provides a powerful framework for studying tissue architecture and gene regulation at single-cell resolution.</p>
