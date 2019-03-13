# elevationAssemblyHawaii

`generateRarefiedMatrices.R`
* Using raw read tables, we generated rarerified read abundances that controlled for the number of individuals in each size category at each site.
    * Each size category at each site were indexed separately
    * For example, if there were 266 individuals between 0 - 2 mm in length, then 2660 reads were rarefied from the raw read pool to get an estimate of abundance while controlling for amplification biases
    * Every read has a zOTU identity (i.e., unique haplotypes)
    * After rarefaction, a master rarefied read abundance matrix is generated (rows: ZOTU IDs, columns: sites)
    * Master site-zOTU matrix was also collapsed into more inclusive levels of hierachy: OTU (zOTUs that are at least 97% similar) and Species (zOTUs that can be clustered)
    * ZOTUs, OTUs and Species that have been identified (using BLAST algorithm) of some orders were omitted from the master matrix (Collembola, Chilopoda, Diplopoda, Blattodea, Isopoda, Mantodea, Phasmatodea) as they are highly likely to be non-native or invasive species in the sample.

`metabarcodingTools.R`
* Script containing functions for rarefaction, collapsing read abundance matrices into more inclusive levels of hierarchy (e.g., zOTU to OTU or OTU to Species)


`analysis.R`
* 