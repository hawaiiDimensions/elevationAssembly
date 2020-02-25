# elevationAssemblyHawaii

`generateRarefiedMatrices.R`
* Using raw read tables, we generated rarerified read abundances that controlled for the number of individuals in each size category at each site.
    * Each size category at each site were indexed separately
    * For example, if there were 266 individuals between 0 - 2 mm in length, then 266 * 7 reads were rarefied from the raw read pool to get an estimate of abundance while controlling for amplification biases
    * 7 reads as that was the maximum number of reads that would work for the sample with the fewest number of reads
    * Every read has a zOTU identity (i.e., unique haplotypes)
    * After rarefaction, a master rarefied read abundance matrix is generated (rows: ZOTU IDs, columns: sites)
    * Master site-zOTU matrix was also collapsed into more inclusive levels of hierachy: OTU (zOTUs that are at least 97% similar) and Species (zOTUs that can be clustered)
    


`arf_rarefied.rds`, `mco_rarefied.rds`
* List of rarefied samples using either the ARF loci or the MCO loci

`combinedOTU_data_rX.rds`
* Number of rarefied reads from both ARF and MCO loci were combined

`prepOTU.R`
* For each combined OTU table, the sum of rarefied reads from each loci for each OTU were summed at the site-level (i.e., across site and size categories) = `master_zOTU_rX.rds`
* ZOTUs, OTUs and Species that have been identified (using BLAST algorithm) of some orders were omitted from the master matrix (Collembola, Chilopoda, Diplopoda, Blattodea, Isopoda, Mantodea, Phasmatodea) as they are highly likely to be non-native or invasive species in the sample. = `zOTU_native_rX.rds` and `OTU_native_rX.rds`
* For the OTUs, individual orders were also saved separately = `OTU_Lepidoptera_r1.rds`

`metabarcodingTools.R`
* Script containing functions for rarefaction, collapsing read abundance matrices into more inclusive levels of hierarchy (e.g., zOTU to OTU or OTU to Species)


`analysis.R`
* 