# Common densoviruses in the human and mammalian virospheres

Densoviruses (subfamily Densovirinae, family Parvoviridae) are small, non-enveloped T=1 icosahedral ssDNA viruses whose linear genomes are typically ~4 to 6 kb and terminate in hairpin structures that act as replication origins (Cotmore et al., 2019; Bergoin and Tijssen, 2000). Genome organization follows the parvoviral pattern with a 5′ nonstructural block encoding NS proteins including the SF3 helicase/endonuclease NS1, and a 3′ capsid block encoding VP proteins, with extensive use of alternative promoters, leaky scanning, and splicing to generate multiple products (Bergoin and Tijssen, 2000; Guérin et al., 2015). At the genus level, densoviruses show meaningful architectural variation: Iteradensoviruses have a monosense transcription strategy and multiple VP size variants (ICTV Iteradensovirus Report, 2023), whereas Ambidensoviruses use an ambisense arrangement with oppositely oriented NS and VP cassettes (Liu et al., 2015). Replication occurs in the nucleus and requires host S-phase functions, producing prominent intranuclear inclusions (Ren et al., 2017). Consistent with that dependence, cellular tropism is skewed toward mitotically active or differentiating tissues in natural hosts, commonly larval, gut, epidermal, or hematopoietic compartments, and infections can be acutely lytic or persistent depending on host and virus lineage (Bergoin and Tijssen, 2000; Ren et al., 2017). ICTV taxonomy recognizes multiple densovirus genera distributed across two deep evolutionary branches within Densovirinae, underscoring their long-standing diversification in invertebrates (ICTV Densovirinae Report, 2023; ICTV Parvoviridae Report, 2023).￼

Evidence for densovirus infection of mammals remains limited, and most detections in mammalian samples are interpreted as environmental or dietary transit from infected invertebrates rather than productive replication. For instance, densovirus sequences recovered from insectivorous bat feces cluster with arthropod densoviruses, consistent with ingestion and passage through the gut (Li et al., 2012). A frequently cited human-linked genome is human CSF-associated densovirus 1 (HuCSFDV1), detected and independently confirmed in cerebrospinal fluid from a patient with anti-NMDA-receptor encephalitis; the genome is highly divergent but falls within iteradensovirus-like lineages (Phan et al., 2016). Additional human associations include metagenomic detection of densovirus-like sequences in plasma from Cameroonian blood donors participating in blood-borne virus surveillance, again without evidence of a sustained clinical syndrome or demonstrated human tropism (Sadeuh-Mba et al., 2023). Densoviruses have also been noted as low-abundance members of the cutaneous DNA virome in longitudinal shotgun metagenomic surveys of healthy skin, including the Oh et al. skin microbiome study, where densovirus reads were cataloged but not a focus of analysis (Oh et al., 2016). Together, these reports suggest that densovirus sequences can be recovered from human-associated samples, but whether they represent rare vertebrate infection, transient carriage, or contamination remains unresolved, and there is still no robust experimental system showing efficient densovirus replication in vertebrate cells in vitro or in vivo (Phan et al., 2016; Cotmore et al., 2019).  ￼

The interpretation of such sporadic detections is changing rapidly because sequence-based interrogation of the Sequence Read Archive (SRA) has become practical at petabase scale. Early efforts such as Serratus used cloud-based ultra-high-throughput alignment of SRA runs to reference sets, effectively bringing targeted, BLAST-like screening to millions of datasets but still requiring substantial compute and being optimized for family-level discovery rather than arbitrary ad hoc queries (Edgar et al., 2021; Serratus Project, 2020). More recent indexing approaches build searchable k-mer or graph representations of SRA content. Logan-search allows users to submit a short query sequence and receive ranked SRA accessions likely to contain matching reads within minutes (LoganSearch GitHub, 2025; Logan Consortium, 2024). MetaGraph similarly indexes public read sets into compressed de Bruijn graph structures enabling rapid presence/absence and abundance queries over petabases of data, with follow-up alignment providing nucleotide-level confirmation (Karasikov et al., 2024; Karasikov et al., 2025). Practically, these frameworks extend BLAST-like functionality to SRA turning what was previously computationally expensive or prohibitive per query into web searches that can be completed in minutes (Karasikov et al., 2024; Logan Consortium, 2024).

Our group has been developing air sampling as a metagenomic surveillance platform for human viruses. Using Thermo Fisher AerosolSense samplers deployed in built environments, we showed that sequence-independent workflows coupled to deep sequencing can recover diverse respiratory and enteric viruses directly from ambient air, including influenza A/C, RSV, seasonal coronaviruses, rhinovirus, SARS-CoV-2, rotavirus, and astroviruses, validating air metagenomics as a practical complement to clinical and wastewater surveillance (Minor et al., 2023). In parallel, air monitoring has been extended to ports of entry: active samplers placed near customs and immigration areas of major U.S. international airports generated nucleic acids that were processed by hybrid-capture enrichment using Illumina RNA Prep with Enrichment and the Respiratory Virus Oligo Panel, enabling metagenomic recovery of dozens of viral taxa and high-quality genomes for SARS-CoV-2, influenza, bocavirus, and seasonal coronaviruses from airport air samples (Air Monitoring in International Airports Consortium, 2025). Together, these studies establish that low-biomass air samples can be rendered sequence-informative by combining high-volume aerosol capture with unbiased or capture-enriched metagenomics, and that the approach scales from local congregate settings to globally connected travel hubs (Minor et al., 2023; Air Monitoring in International Airports Consortium, 2025).  ￼

During our continued analysis of air metagenomic datasets generated in these efforts, we repeatedly observed densovirus contigs closely related to the previously described HuCSFDV1, appearing across multiple independent air samples and settings. This unexpected recurrence in air samples provided a rationale to revisit the question of densoviruses in mammals using orthogonal data streams. Specifically, it motivated three linked follow-ups: first, systematic mining of public SRA datasets for reads related to HuCSFDV1-like sequences using petabase-scale search tools; second, re-analysis of the longitudinal skin metagenomes in Oh et al. to quantify and contextualize densovirus signal in the cutaneous virome; and third, targeted screening of additional air samples collected through our ongoing surveillance network to map prevalence, diversity, and geographic distribution of densovirus reads in air.

## Results

### HuCSFDV1 is commonly found in air samples from congregate spaces

We collected air samples from congregate indoor spaces in K-12 schools, essentially as previously described. In most sites, a ThermoFisher AerosolSense is run continuously for at least one day. Particulates are eluted in phosphate-buffered saline with 0.1% Tween. Total nucleic acids are extracted from the eluate on a Promega Maxwell instrument. Nucleic acids are enriched for viruses using the Illumina VSP2 protocol according to the manufacturer's recommendation. The purpose of target enrichment is for monitoring of conventional human pathogens. Densoviruses are not included in the VSP probes and are therefore not specifically enriched. Sequences were obtained on an Illumina NovaSeqX. Three datasets with [102 total samples](air-samples/vsp/vsp-sample-metadata.md) were sequenced in three independent VSP2 experiments.

We developed a [Snakemake workflow](data_processing/read_mapping) aligns reads to the HuCSFDV1 reference [NC_076998](https://www.ncbi.nlm.nih.gov/nuccore/NC_076998.1) using minimap2 (v2.28-r1209) with short-read presets, and the resulting alignments are processed through samtools (v1.20) for coordinate sorting, mate-pair fixing, and PCR duplicate removal, yielding indexed BAM files suitable for variant calling and coverage analysis.

76 of the 102 samples from K-12 schools [had at least one read mapped to HuCSFDV1](https://dholab.github.io/common-densoviruses/dane_county_air/). While spurious read mapping could be explained by non-specific mapping, 14 of the samples had at least 100 reads mapped to HuCSFDV1. The sample with the strongest read support for HuCSFDV1, termed `high-school-001-20250901-pooled`, was collected between 2025-09-04 and 2025-09-08. 2,355 reads from this sample mapped to HuCSFDV1. The majority consensus sequence from this sample is 99% nucleotide identical to HuCSFDV1, with 48 majority consensus nucleotide substitutions.

Fortuitously, two sites (high-school-001 and elementary-school-004) each yielded more than 100 mapped reads in independent early, mid, and late September 2025 air samples. For each September collection, we generated consensus sequences from indexed BAM files with `samtools consensus` (v1.20) using FASTA output and retaining ambiguous sites as N to avoid biasing downstream distance calculations, and [aligned these sequences](air-samples/vsp/align/consensus_plus_ref.aln.fasta) together with the NC_076998 reference using MAFFT (v7.526) with default nucleotide parameters. Phylogenies were inferred with IQ-TREE (v3.0.1) under a GTR+G model with 1,000 ultrafast bootstrap replicates, rooting on NC_076998 to provide an external outgroup. Across trees, the three elementary-school-004 consensuses [formed a monophyletic cluster distinct from the three high-school-001 consensuses](air-samples/vsp/tree/consensus_plus_ref.treefile), which similarly grouped together, consistent with independent circulation of related but distinct viruses in the two schools rather than repeated sampling of a shared source.

We recently participated in another air sampling project in [international airports](https://www.medrxiv.org/content/10.1101/2025.09.22.25336185v3). Metagenomic sequencing with VSP enrichment was also performed on a subset of these samples (Bioproject: PRJNA989177), resulting in 159 datasets. Reads from these datasets were downloaded and mapped to the HuCSFDV1 reference as described above. Seven samples had a [single paired-end read](air-samples/tgs/consensus) that mapped to HuCSFDV1; a BLASTN search against core-nt (as of 23 November 2025) unambiguously HuCSFDV1 as the most significant alignment for all seven. Unlike schools where air samplers are located primarily in cafeterias where students and staff dwell, the air samplers in airports were located in international arrival corridors. Arriving passengers would be in the vicinity of the air samplers briefly, possibly accounting for the less robust detection of HuCSFDV1 reads. Nevertheless, this demonstrates that reads corresponding to this virus were found in multiple air samples in two international airports, Dulles and San Fransisco, separated by thousands of miles.

### HuCSFDV1 is detected in mammalian datasets in NCBI SRA

The HuCSFDV1 reference [NC_076998](https://www.ncbi.nlm.nih.gov/nuccore/NC_076998.1) was split into 1kb segments. Each segment was input into [Logan Search](https://logan-search.org/dashboard) using the default treshold of 0.5 and the 'All' reference group. [The results](logan/logan-output) of Logan Search for each segment were combined to create the set of matching SRA accessions. [This table](logan/sra-metadata.md) shows all of the matching datasets. If the source of HuCSFDV1 is a reagent contaminant we would expect to see datasets with matching kmers throughout SRA, largely irrespective of sample origin. Conversely, if HuCSFDV1 is an invertebrate virus like other densoviruses, we would expect the most matches in invertebrate SRA datasets. This is not what we observed. Instead, of the 126 SRA datasets with matches sequences, 56 are from mice. The next most common data sources are air and wastewater genomes, with 16 and 10 datasets, respectively. 5 human datasets also have matches to HuCSFDV1, including XX datasets that were deposited as part of the original description of HuCSFDV1 in 2016.

Reads from all of these SRA datasets were downloaded with fasterq-dump, mapped to the HuCSFDV1 reference and deduplicated as described above, yielding indexed BAM files. 

### HuCSFDV1 is found in human skin microbiome

In a 2016 manuscript, metagenomic sequencing was performed 594 skin samples from 12 healthy individuals. A subset of these individuals had reads that were categorized as Acheta domestica (i.e., house cricket) densovirus. Sequencing reads from the manuscript are available in [NCBI BioProject 46333](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=46333). Unfortunately, there are 4,619 SRA datasets in this BioProject, along with 5,894 BioSamples. The umbrella BioProject does not list this publication in its references, however, exploration of the BioProject identified [1,120 metagenomic datasets](oh-et-al-skin-microbiome/sra-experiment-records.md) with identifiers that matched the ones used in the manuscript (SRA Links for BioProject (Select 46333) AND "human skin metagenome"[orgn:__txid539655]).

Additionally, [Supplemental Table S2 from Oh et al.](https://www.cell.com/cms/10.1016/j.cell.2016.04.008/attachment/db43522c-ed59-47f2-a378-469abee3924b/mmc3.xlsx) describes 585 samples with taxonomic classifications, including 56 samples with densovirus reads. We reformatted this into [a table](oh-et-al-skin-microbiome/samples-with-densovirus-status.md) that describes the 585 samples, the corresponding data in NCBI, and densovirus status.

that accepts BioSample accessions as input, automatically retrieves all associated SRA run accessions via Entrez queries, downloads reads using fasterq-dump (v3.1.1), and concatenates multi-run samples into unified FASTQ files. Reads are aligned 



## Authors and affiliations

Isla Emmen<sup>1</sup>, Andrew Lail<sup>1</sup>, Nancy A. Wilson<sup>1</sup>, Will Gardner<sup>1</sup>, Nicholas R. Minor<sup>1</sup>, Marc C. Johnson<sup>2,3</sup>, Shelby L. O'Connor<sup>1,4</sup>, David H. O'Connor<sup>1,4</sup>

<sup>1</sup> Department of Pathology and Laboratory Medicine, University of Wisconsin School of Medicine and Public Health, Madison, WI 53711, USA
<sup>2</sup> Department of Molecular Microbiology and Immunology, University of Missouri School of Medicine, Columbia, MO 65211, USA
<sup>3</sup> Christopher S. Bond Life Sciences Center, University of Missouri, Columbia, MO 65211, USA
<sup>4</sup> Wisconsin National Primate Research Center, University of Wisconsin-Madison, Madison, WI 53715, USA

## Acknowledgements

This observation would not have been possible without data submitted to NCBI SRA by Oh and colleagues in Bioproject 46333, as well as the other SRA datasets where we detected densovirus reads. A full list of these datasets and their submitters is in XX.

Searching SRA datasets for densovirus reads was made possible by [Logan Search](https://logan-search.org) ([Github](https://github.com/IndexThePlanet/LoganSearch)), and the work of the Logan Team.

This work is funded by Inkfish.

Computing was performed on the UW-Madison Center for High Throughput Computing cluster. Claude Code, Claude Sonnet v4.5, Claude Opus v4.5, and ChatGPT v5.1 were used for research and text drafting.

## References

Bergoin M, Tijssen P. Biological and molecular properties of densoviruses and their use in protein expression and biological control. In: Miller LK, Ball LA, editors. The Insect Viruses. Springer; 2000. p. 141–169.

Cotmore SF, Agbandje-McKenna M, Canuti M, et al. ICTV Virus Taxonomy Profile: Parvoviridae. J Gen Virol. 2019;100:367-368. doi:10.1099/jgv.0.001212.

Edgar RC, Taylor J, Lin V, et al. Petabase-scale sequence alignment catalyses viral discovery. bioRxiv. 2021. doi:10.1101/2020.08.07.241729.

Guérin F, Lapointe R, Sénéchal F, et al. Gene expression of five different iteradensoviruses. J Virol. 2015;89:1955-1968. doi:10.1128/JVI.01719-14.

Hewson I, Button JB, Gudenkauf BM, et al. Densovirus associated with sea-star wasting disease and mass mortality. Proc Natl Acad Sci U S A. 2014;111:17278-17283. doi:10.1073/pnas.1416625111.

Hewson I, Bistolas KSI, Quijano Cardé EM, et al. A highly prevalent and pervasive densovirus discovered among sea stars from the North American Atlantic coast. Appl Environ Microbiol. 2019;85:e02723-19. doi:10.1128/AEM.02723-19.

ICTV. Subfamily Densovirinae. In: ICTV Report: Parvoviridae. 2023 release (MSL #39).

Karasikov M, et al. Indexing and searching petabase-scale nucleotide resources. Nat Methods. 2024;21:xxx-xxx. doi:10.1038/s41592-024-02280-z.

Karasikov M, et al. Efficient and accurate search in petabase-scale sequence repositories (MetaGraph). Nature. 2025;xxx:xxx-xxx. doi:10.1038/s41586-025-09603-w.

Ledermann JP. Densonucleosis viruses (densoviruses) for mosquito and pathogen control. Curr Opin Insect Sci. 2018;28:57-63. doi:10.1016/j.cois.2018.05.006.

Li L, Victoria JG, Wang C, et al. Bat guano virome: predominantly insect viruses and a densovirus genome. J Virol. 2012;86:4620-4627. doi:10.1128/JVI.06671-11.

Liu W, Zhang Q, Yang B, et al. First complete genome of an ambidensovirus from crayfish Cherax quadricarinatus and evidence for ambisense transcription strategy. Virus Res. 2015;208:55-62. doi:10.1016/j.virusres.2015.05.012.

Logan Consortium. Logan: planetary-scale genome assembly surveys life’s diversity. bioRxiv. 2024. doi:10.1101/2024.07.30.605881.

LoganSearch GitHub. Logan-search: k-mer search engine for all SRA public accessions. 2025. https://github.com/IndexThePlanet/LoganSearch.

Meynadier G, Vago C, Plantevin G. Isolement d’un virus de type parvovirus chez les larves de la grande fausse-teigne Galleria mellonella et description de la densonucleose. Rev Zool Agric Appl. 1964;63:207-209.

Oh J, Byrd AL, Deming C, et al. Temporal stability of the human skin microbiome. Cell. 2016;165:854-866. doi:10.1016/j.cell.2016.04.008.

Phan TG, Messacar K, Dominguez SR, da Costa AC, Deng X, Delwart E. A new densovirus in cerebrospinal fluid from a case of anti-NMDA-receptor encephalitis. Arch Virol. 2016;161:3231-3235. doi:10.1007/s00705-016-3002-9.

Sadeuh-Mba SA, et al. Metagenomic detection of divergent insect- and bat-associated viruses in plasma from two African individuals enrolled in blood-borne surveillance. Viruses. 2023;15:1022. doi:10.3390/v15041022.

Serratus Project. Ultra-deep search for novel viruses. 2020–2025. https://serratus.io/.
