#' proteins
#'
#' A dataset containing human proteome information from Uniprot.
#'
#' @usage data(proteins)
#'
#' @docType data
#'
#' @format A data frame with 20,430 rows and 8 variables:
#' \describe{
#'   \item{uniprot_id}{Official uniprot entry id}
#'   \item{gene_name}{Gene name}
#'   \item{gene_name_alt}{Alternative gene names associated with entry}
#'   \item{protein_name}{Protein name}
#'   \item{protein_name_alt}{Alternative protein names associated with entry}
#'   \item{sequence}{Primary amino acid sequence}
#'   \item{length}{Number of amino acids}
#'   \item{mass}{Molecular weight of protein}
#'   }
#'
#' @keywords datasets
#'
#' @source \url{https://www.uniprot.org/uniprot/?query=*&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+reviewed%3Ayes}
"proteins"

#' Human chromosome summary
#'
#' A dataset containing summary information of the human chromosomes
#'
#' @usage data(chromosome)
#'
#' @docType data
#'
#' @format A data frame with 24 rows and 14 variables:
#' \describe{
#'   \item{id}{Identified for each of the 23 chromosome pairs in cell nuclei, including the X and Y chromosome; mitochondrial DNA is excluded.}
#'   \item{length_mm}{Chromosome length calculated by multiplying the number of base pairs by 0.34 nanometers, the distance between base pairs in the DNA double helix.}
#'   \item{basepairs}{Number of base pairs, based on a reference genome that does not represent the sequence of any specific individual (Data source: Ensembl genome browser release 87, December 2016 for most values).}
#'   \item{variations}{Number of variations of unique DNA sequence differences that have been identified in the individual human genome sequences analyzed by Ensembl as of December 2016.}
#'   \item{protein_codinggenes}{Number of enes encoding proteins is based on the number of initial precursor mRNA transcripts, and does not include products of alternative pre-mRNA splicing, or modifications to protein structure that occur after translation.}
#'   \item{pseudo_genes}{Number of non-functional remnants of protein-coding genes}
#'   \item{totallongnc_rna}{Number of long non-coding RNAs -- RNA molecules longer than 200 bases that do not have protein-coding potential. These include: ribosomal RNAs, or rRNAs (the RNA components of ribosomes), and a variety of other long RNAs that are involved in regulation of gene expression, epigenetic modifications of DNA nucleotides and histone proteins, and regulation of the activity of protein-coding genes.}
#'   \item{totalsmallnc_rna}{Number of small non-coding RNAs -- RNAs of as many as 200 bases that do not have protein-coding potential. These include: microRNAs, or miRNAs (post-transcriptional regulators of gene expression), small nuclear RNAs, or snRNAs (the RNA components of spliceosomes), and small nucleolar RNAs, or snoRNA (involved in guiding chemical modifications to other RNA molecules).}
#'   \item{mi_rna}{Number of microRNA (abbreviated miRNA) -- small non-coding RNA molecule (containing about 22 nucleotides) found in plants, animals and some viruses, that functions in RNA silencing and post-transcriptional regulation of gene expression.}
#'   \item{r_rna}{Number of ribosomal ribonucleic acid (rRNA) -- the RNA component of the ribosome, which is essential for protein synthesis in all living organisms. rRNA is the predominant RNA in most cells, composing around 80% of cellular RNA. Ribosomes are approximately 60% rRNA and 40% protein by weight. A ribosome contains two subunits, the large ribosomal subunit (LSU) and small ribosomal subunit (SSU).}
#'   \item{sn_rna}{Number of small nuclear RNA (snRNA) -- a class of small RNA molecules that are found within the splicing speckles and Cajal bodies of the cell nucleus in eukaryotic cells. The length of an average snRNA is approximately 150 nucleotides.}
#'   \item{sno_rna}{Number of small nucleolar RNAs (snoRNAs) -- a class of small RNA molecules that primarily guide chemical modifications of other RNAs, mainly ribosomal RNAs, transfer RNAs and small nuclear RNAs. There are two main classes of snoRNA, the C/D box snoRNAs, which are associated with methylation, and the H/ACA box snoRNAs, which are associated with pseudouridylation.}
#'   \item{miscnc_rna}{Number of miscellaneous RNA, a general term for a series of miscellaneous small RNA. It serves a variety of functions, including some enzyme-like catalysis and processing RNA after it is formed.}
#'   \item{centromereposition_mbp}{Numeric indicating the position of the region that joins the two sister chromatids, or each half of the chromosome.}
#' }
#'
#' @keywords datasets
#'
#' @source \url{https://en.wikipedia.org/wiki/Human_genome}
"chromosome"

#' gene_location
#'
#' A dataset containing information of human genes and their location across the human genome
#'
#' @usage data(gene_location)
#'
#' @docType data
#'
#' @format A data frame with 20388 rows and 14 variables:
#' \describe{
#'   \item{ensembl_transcript_id}{Ensembl transcript id}
#'   \item{ensembl_gene_id}{Ensembl gene id}
#'   \item{chromosome_name}{Chromosome name}
#'   \item{transcript_start}{Transcript start site}
#'   \item{transcript_end}{Transcript end site}
#'   \item{strand}{Strand}
#'   \item{band}{Band}
#'   \item{approved_symbol}{Approved gene symbol}
#'   \item{percentage_gene_gc_content}{Percent GC content}
#'   \item{gene_name}{Gene name}
#'   \item{description}{Gene description}
#'   \item{cds_length}{CDS length}
#' }
#'
#' @keywords datasets
#'
#' @source \url{http://uswest.ensembl.org/index.html}
"gene_location"
