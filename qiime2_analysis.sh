#!/bin/bash

# ============================================================
# Gut Microbial Dysbiosis and Elevated Mucosal Inflammation
# in Severe Acute Malnutrition (SAM)
#
# QIIME 2 analysis workflow
# QIIME 2 environment: qiime2-amplicon-2024.10
# ============================================================

qiime demux summarize \
  --i-data demux-trimmed.qza \
  --o-visualization demux-trimmed.qzv

  qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-trimmed.qza \
  --p-n-threads 6 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 227 \
  --p-trunc-len-r 224 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats stats.qza

mkdir -p dada2_output
mv table.qza rep-seqs.qza stats.qza dada2_output/

qiime feature-table summarize \
  --i-table dada2_output/table.qza \
  --o-visualization table.qzv

qiime feature-table tabulate-seqs \
  --i-data dada2_output/rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file dada2_output/stats.qza \
  --o-visualization stats.qzv

qiime feature-classifier classify-sklearn \
  --i-classifier '/home/raj/Downloads/silva-138-99-nb-classifier.qza' \
  --i-reads dada2_output/rep-seqs.qza \
  --o-classification dada2_output/taxonomy.qza

qiime metadata tabulate \
  --m-input-file dada2_output/taxonomy.qza \
  --o-visualization dada2_output/taxonomy.qzv

  qiime taxa filter-table \
  --i-table dada2_output/table.qza \
  --i-taxonomy dada2_output/taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table dada2_output/table-no-mito-chloro.qza

qiime taxa filter-seqs \
  --i-sequences dada2_output/rep-seqs.qza \
  --i-taxonomy dada2_output/taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences rep-seqs-no-mito-chloro.qza

  qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences dada2_output/rep-seqs-no-mito-chloro.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table dada2_output/table-no-mito-chloro.qza \
  --p-sampling-depth 40000 \
  --m-metadata-file dada2_output/metadata.tsv \
  --output-dir core-metrics-results

  # Shannon diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file dada2_output/metadata.tsv \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

# Faith's Phylogenetic Diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file dada2_output/metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-significance.qzv

# Observed features
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file dada2_output/metadata.tsv \
  --o-visualization core-metrics-results/observed-features-significance.qzv

  # Bray-Curtis
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file dada2_output/metadata.tsv \
  --m-metadata-column Group \
  --p-method permanova \
  --o-visualization core-metrics-results/bray-curtis-permanova.qzv \
  --p-pairwise

# Unweighted UniFrac
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file dada2_output/metadata.tsv \
  --m-metadata-column Group \
  --p-method permanova \
  --o-visualization core-metrics-results/unweighted-unifrac-permanova.qzv \
  --p-pairwise

# Weighted UniFrac
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file dada2_output/metadata.tsv \
  --m-metadata-column Group \
  --p-method permanova \
  --o-visualization core-metrics-results/weighted-unifrac-permanova.qzv \
  --p-pairwise

  # Phylum level
qiime taxa collapse \
  --i-table dada2_output/table-no-mito-chloro.qza \
  --i-taxonomy dada2_output/taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table phylum-table.qza

# Family level
qiime taxa collapse \
  --i-table dada2_output/table-no-mito-chloro.qza \
  --i-taxonomy dada2_output/taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table family-table.qza

# Genus level
qiime taxa collapse \
  --i-table dada2_output/table-no-mito-chloro.qza \
  --i-taxonomy dada2_output/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table genus-table.qza

mkdir -p exports

# ASV table
qiime tools export \
  --input-path dada2_output/table-no-mito-chloro.qza \
  --output-path exports

mv exports/feature-table.biom exports/asv-table.biom

biom convert \
  -i exports/asv-table.biom \
  -o exports/asv-table.tsv \
  --to-tsv


# Phylum table
qiime tools export \
  --input-path phylum-table.qza \
  --output-path exports

mv exports/feature-table.biom exports/phylum-table.biom

biom convert \
  -i exports/phylum-table.biom \
  -o exports/phylum-table.tsv \
  --to-tsv


# Family table
qiime tools export \
  --input-path family-table.qza \
  --output-path exports

mv exports/feature-table.biom exports/family-table.biom

biom convert \
  -i exports/family-table.biom \
  -o exports/family-table.tsv \
  --to-tsv


# Genus table
qiime tools export \
  --input-path genus-table.qza \
  --output-path exports

mv exports/feature-table.biom exports/genus-table.biom

biom convert \
  -i exports/genus-table.biom \
  -o exports/genus-table.tsv \
  --to-tsv
