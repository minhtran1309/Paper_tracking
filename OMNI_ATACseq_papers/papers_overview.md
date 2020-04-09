# Analysis of data for 2 initial papers

There are two papers for that project:
1. [Overview](#overview)
  1. [Ageing "roadmap"](#project1)
  2. [Predictive framework for cell type/state interconversion](#project2)
2. [Specific of Project 2](#spec_proj_2)
3. [References](#references)
## <a id="overview">1. Overview</a>
-----
### <a id="project1">a. Ageing "roadmap"</a>
---
1. Characterise TF network changes in Bulk sapmles
2. Characterise TF network changes in sc samples. **For hematopietic system**
3. Verify reallocation of key TKs with CutnRun of CutnTag for 2-3 key cell types (possibly integrate with DNA methylation data). 
4. Perturb key TFs in cell type (e.q. T-cells or HSCs) and assess consequences at the end of the year (scRNAseq + scATACseq)

### <a id="project2">b. Predictive framework for cell type/state interconversion</a>
---
1. Categorise TFs across the data set into activate and repressor (Diff TF)
2. Characterise TF activity within and across cell states (DiffTF, Hint seq, Differential peak calling + Motive analysis). Indentify TFs with strongest/pioneer activity
3. Determine "Reguatory" peaks (peaks that correlate with gene expression +/- 500kb of TSS of a gene-use entire dataset).
4. Perform "Mogrify-like" predictions based on information about (I) regressor/Activator function and (II) use of genes' regulatory peaks (II) for motive prediction (rather than just 500bps surrounding TSS like Mogiry).
5. Adjust TF prediction/scores (state A to B) with TF activity results between A and B (as per III) - ideally priotize TFs with very strong or pioneer activity.
---
## <a id="spec_proj_2">2. Specific about project 2</a>
Priotize project 2
1. For the first goals of that project: refer to [Berest et al., 2019](#Berest_et_al_2019), Cell Reports
  * In short: Adapt DiffTF to identify TFs that are activators and repressors by correlating their expression lvl with accessibility at their TFs binding sites across the entire dataset.
  * In detail: This is the original work of DiffTF
    * What is DiffTF about? -> a tool to estimate differential TF activity and classify TFs into activators or repressors.
    * What DiffTF needs as input? -> active chromatin data (accessibity/ChIP-seq) and **integrates** with RNA-seg for classification.
    * **Summary** -> DiffTF, in basic mode, calculates differential TF activity, combines genome-wide chromatin accessibity/activity with putative TF binding sites that. In classification mode, TF binding sites are integrated with RNA-seq. It's applied to compuate mutated and unmutated chronic lymphocytic leukemia patients and two hematopoietic progenitor cell types.
    * How do we define classes of TFs? -> Since most TFs have been reported to act as both activator and repressor, depending on the study, most likely that there is no definite term for classifying activator and repressor. Therefore, the classification framework is based on the following assumption: (1) increasing the abundance of an **activator TF** will **increase the average accessibility** at the regulatory elements controlled by the TF which will **lead to upregulation of its target genes**. (2) increasing the abundance of a **repressor TF** will **decrease the average accessibility** at the regulatory elements controlled by the TF, which will lead to downregulation of its target genes. The second assumption about repressor may not be straighforward. 
    * There are two case studies: **Case I**, Differential TF activity in a Heterogeneous ATAC-Seq Dataset in CLL. The experiment carried out on ATAC-seq data for 56 CLL patients (88 samples) that are statified by the mutation (originally there are 34 U-CLL and 50 M-CLL and 4 unclassified), after data processing and QC, (25 U-CLL and 27 M-CLL remained). False discovery rate <10%, which means quite good. **Case II**, Applying diffTF to Small-scale multiomics dataset. The study design is to assess the differences between two condtions with few biological replicates per condition. They used eight ATAC-seq and RNA-seq profiles of multipotent progenitor cells, an early hematopoiecteic progenitor population capable of supporting multilineage blood production
    * How does it work in detail? -> calculate differential TF activity between two or more conditions for each TF by comparing the distribution of fold-change differences across all binding sites of a TF to all binding sites from all other TFs. The algorithm is split into 7 step as described in the following:
         1. Generate consenus peak set: analyze read counts peaks, which requires a consensus peak set across input samples.
For the later, consensus peaks are generated with the function in Bioconductor package. Retain only peaks from genuine autosomes, thereby filtering sex chromosomes, non-assembled contigs as well as alternative haplotypes. The peak set is finally sorted by coordinate to speed up subsequent computations. 
         2. Scanning of TF binding sites, for each TF of interest, diffTF needs a set of TF binding sites. They used the a database called HOCOMOCO, which provides TF binding models.
Finally, they sort the TFBS for each TF by coordinate.
         3. Differential analysis for the consensus peakset, to calculate the fold change between two condition across each peak, they **first** obtain the counts for the consensus peakset for each sample using featureCounts from the Subread package with the options. **Then** they employ DESeq2 for count normalization using a cycling loess approach as implemented for normOffsets function of csaw from Bioconductor. Unrestricted design formulas can be incorporated into diffTF, making it a flexible approach. DiffTF can also incorporate design formulas for which the predictor varibal is continuously-valued rather than binary to analyze time-course data. 
        4. Signal extraction for each TFBS: they filter for TFBS that overlab with the consensus peak set only while allowing for multiple TFBS per peak (using the program bedtools intersect with options). Each TFBS is then extended by 100bp (adjustable) in both directions followed by extraction of read counts for each sample. 
        5. Calculate of accessibility fold-change for each TFBS: this step is to avoid biases and dependencies based on TFBS clustering within peaks, they select the TFBS per TF per peak with highest average read count accross all samples. The count normalization, we use the consensus peak-derived normalization factors from step 3. The result of this step is log2 fold-change value for each selected TFBS per TF per peak, and in addition, are the various diagnostic plots for each TF. 
        6. GC binning and calculation of differential TF activity values: GC binning to reduce the biological biases based on differential effects depending on GC content of the local environment of a TF, using bedtools. Next is the estimation of differntial TF activity, for each bin containing more than 20 TFBS for a given TF. We compare the TF-specific distribution of log2 fold-changes against a background of log2 fold-changes of all TFBS from all other TFs of the same GC bin. Finally, for each TF, the differential TF activity is calculated as the weighted arithmetic mean across all mean difference values for the bins with sufficient data, weighting each value from each bij by the fraction of TFBS it contain in that TF such that all weight sum to 1. 
       7. Estimation of significance for differential activity for each TF: to assess their statistical significance, we emplow a permutation approach to derive empirical p values. They rerun steps 3-6 for a total of 1000 times (adjustable) with permuted input data and then calculate an empirical two-sided p values per TF by comparing the real value with the distribution from the permuations and calculating the proportion of sampled permutation for which the absolute differential TF activity is larger. Estimation of the variance is calculated by linear combinations of weighter mean T-scores for each TF. Subsequently, to obtain a p value for each TF, they centralize the distribution of the weighted T scores across all TF by subtracting its mean. They calculate p values out of z-scores based on the TF-specific variance calulated above.  
  * What about classification mode? -> They tried to on literature-mining, using TRRUST -> **Failed, since TFs were classified almost equally often as activator and repressor**. They propose an assumption that increasing the level of an activating TF increases chromatin accessibility at its target sites while increasing the level of a repressing TF decrease it. So based on that assumption, they normalized RNA-seq count data (quantile normalized to minimize the effect of outlier values), then they calculate the pearson correlation coefficients between the expression level of each TF and the ATAC-seq signal of each putative TFBS across all individual 
  * A few notes on ATAC-seq processing: Snakemake pipeline that starts with raw fastq, goes through multiple steps for quality control, adaptor trimming, alignment, as well as general and ATAC-seq specific post-alignment filtering, preprocessing. FastQC is used to assess the seq quality. Trimmomatic for removing the Nextera Transposese agent, alligned with Bowtie2, followed by various cleaning steps (picard tool CleanSam, FixMateInformation. Various filtering steps include: removing mitochondrial read and reads from non-assembled contifs or alternative haplotypes, low quality mapping, remove duplicate read. Lastly, a GC bias diagnosis and correction using deepTools and Benjamin's method is run for each sapmle. 

2. Approaches like DiffTF ("basic mode") and Hint-Seq (Footprinting), refer to [Li et al](#Li_et_al_2019).
  * In short: To identify the TFs with strongest/pioneer activity
  * In detail: 
    * What is Hint-seq about? -> Hint-seq is a footprinting methods considering ATAC-seq protocol artifact, it helps with predicting the transcription factor binding site using footprints by answering the floowing questions:
    * Why is Hint-seq useful? -> ATAC-seq allows the detection of open chromatin by identifying genomic intervals with many reads. However, the presence of TFs bound to DNA prevent the enzyme from cleavage in an otherwise nuclesome-free region, which leaves small regions, **aka footprints, where read coverage suddenly drops within peak region of high coverage**.
    * What Hint-seq needs as input? -> a give open chromatin library and a reference genome sequence G with length N. 
    * Steps summary: First, genomic cleavage signal are generated from raw seq libraries after filtering reads by fragment size, correction of cleavage bias and signal normalization. Second, cleavage are given as input to a HMM, which segment the signal and finds the local of footprints.
    * How does it work in detail?
      1. Cleavage event counting and correction of seq bias: Given ATAC-seq read aligned with start position i, HINT considers the postion i+4 (because the middle of Tn5 claevage event is the fifth base after the fragment start) as a cleavage event for forward reads and i-5 for reverse reads. Next, the correction of the cleavage event profiles by sequence-specific cleavate bias considering the word w[i] with size k around genomic position i. 
      2. k-mer-based estimation: the most common approach for bias estimation is to use freq of k-mers to estimate the probability p(w|obs). 
      3. PWM-based estimation: and alternative approach to calculate bias estimation. p(w|obs).
      4. PDM-based estimation: yet another approach to calucate bias estimation using a given multiset W, to estimate p(wj).
      5. HMM training and decoding: they take the previously describe multivariate cleavage signals X as input for HMM model. 
3. Not so clear, refer to **2 papers**: Ucar et al., 2017 GEM. Moskowits et al., 2017 Science immunology. 

4. Firstly, from RNA-seq produce DE genes for clusters (assuming Cell type A, Cell type B). Secondly, from ATAC-seq produce differential accessible regulatory peaks (based on regulatorypeaks for expressed transcript whole data set), so that we can have list of genes from  Cell type A/B which are more accessible. Finally, project the result from two analyses and perform Mogrify like prediction. 

5.  

## <a id="references">3. References</a>
<a id="Berest_et_al_2019">[Berest et al.]</a> Ivan Berest, Quantification of Differential Transcription Factor Activity and Multiomics-Based Classification into Activators and Repressors: DiffTF, Cell Reports, 2019.
<a id="Li_et_al_2019">[Li et al.]</a> Zhijian Li, Identification of transcription factor binding sites using ATAC-seq: Hint-Seq, Genome Biology, 2019 
