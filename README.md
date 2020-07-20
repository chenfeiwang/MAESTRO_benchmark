# Benchmark codes used in MAESTRO paper
We provided the codes used in evaluating ```scATAC-seq clustering```, ```automatic cell-type annotation```, ```integration between scRNA-seq and scATAC-seq using different methods```, and ```integration evaluation using different peak-RP model``` in this repository. While the dataset for evaluation is a big large, users can download the data at [here](http://cistrome.org/~chenfei/MAESTRO/MAESTRO_benchmark.tar.gz), and put the data in the same directory with the benchmark codes.

MAESRO paper can be found at [citation](citation)

### Clustering benchmark
The data and R codes can be found at ```clustering_benchmark``` directory. Including three types of comparisons.
1) ```clustering_benchmark_methods_simulated_data```
ScATAC-seq clustering using scABC, cisTopic, snapATAC, LSI-based(MAESTRO) on simulated dataset from 10 different bulk ATAC-seq dataset. We sampled the reads at 1000, 2500, 5000, and 10000. The clustering performance was evaluated using [Nomalized Mututal Information](https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html) between cluster labels and original sample type 

2) ```clustering_benchmark_methods_public_data``` 
ScATAC-seq clustering using scABC, cisTopic, snapATAC, LSI-based(MAESTRO) on published microfluidic cellline, HSC and 10X Genomics PBMC dataset. For cellline and HSC dataset, clustering accuracy was evaluated using NMI between cluster labels and original cell labels. For PBMC, clustering accuracy was evaluated using [RAGI](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1854-5).

3) ```clustering_benchmark_peaks```
ScATAC-seq clustering using different peak sets, single-cell peak only, ENCODE CCRE only or single-cell peak combined with ENCODE CCRE.For cellline and HSC dataset, clustering accuracy was evaluated using NMI between cluster labels and original cell labels. 

### Cell-type annotation benchmark
The data can R codes be found at ```annotation_benchmark``` directory.
We benchmarked the performance of automatic celltype annotation using SCINA, Garnett and MAESTRO on sorted PBMC dataset from [Zheng et,al](https://www.nature.com/articles/ncomms14049). The cell-types were annotated using LM22 or a simple gene signature from Garnett. The performance was evaluated using median F1-score between true labels and annotated cell-type labels. Codes for benchmarking were adopted from [scRNAseq_Benchmark](https://github.com/tabdelaal/scRNAseq_Benchmark).

### Integration between scRNA-seq and scATAC-seq using different methods
The data can R codes be found at ```integration_benchmark``` directory.
We benchmarked the integration performance using MAESTRO peak-RP model, Seurat gene activity score, snapATAC-seq genebody score, and cicero score. The performance was evaluated by the distribution of cell-type label prediction score, and spearman correlation between scRNA-seq and scATAC-seq gene activity score on both top 2000 highly variable genes and all genes.

### Integration between scRNA-seq and scATAC-seq using different peak-RP model
The data can R codes be found at ```integration_RPmodel_benchmark``` directory.
We benchmarked the integration performance of different MAESTRO peak-RP models. The performance was evaluated by the distribution of cell-type label prediction score, and spearman correlation between scRNA-seq and scATAC-seq gene activity score on top 2000 highly variable genes.