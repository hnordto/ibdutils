
snp2vcf <- function(SNPdata,
                     filename = NULL) {

  vcf <- data.frame(CHROM = SNPdata$CHROM,
                    POS = SNPdata$MB*10^6,
                    ID = SNPdata$MARKER,
                    REF = SNPdata$A1,
                    ALT = SNPdata$A2,
                    QUAL = rep(".",nrow(SNPdata)),
                    FILTER = rep("PASS",nrow(SNPdata)),
                    INFO = paste0("AF=",1-SNPdata$FREQ1),
                    FORMAT = rep("GT",nrow(SNPdata)),
                    ID1 = apply(SNPdata, 1, genotypeToNumeric, genotypeCol = 8),
                    ID2 = apply(SNPdata, 1, genotypeToNumeric, genotypeCol = 9))
  header <- c(
    "##fileformat=VCFv4.2",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
    "##FORMAT<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tID1\tID2"
  )

  if (!is.null(filename)) {
    writeLines(header, con = filename)
    write.table(vcf, file = filename, append = TRUE, quote = FALSE,
                sep ="\t", row.names = FALSE, col.names = FALSE)
  }

  return (vcf)
}

# Only supporting unphased data as of now

genotypeToNumeric <- function(data, genotypeCol) {

  genotype_str <- data[genotypeCol]

  dict = list()
  dict[[data[5]]] = 0
  dict[[data[6]]] = 1

  genotype_str_vec <- unlist(strsplit(genotype_str,"/"))

  genotype_num <- paste0(dict[[genotype_str_vec[1]]],
                         "/",
                         dict[[genotype_str_vec[2]]])
  return(genotype_num)



}


