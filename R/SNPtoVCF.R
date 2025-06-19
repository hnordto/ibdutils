
SNPtoVCF <- function(SNPdata,
                     vcf_file) {
  data <- get(load(SNPdata))

  vcf <- data.frame(CHROM = data$CHROM,
                    POS = data$MB,
                    ID = data$MARKER,
                    REF = data$A1,
                    ALT = data$A2,
                    QUAL = rep(".",nrow(data)),
                    FILTER = rep("PASS",nrow(data)),
                    INFO = paste0("AF=",1-data$FREQ1),
                    FORMAT = rep("GT",nrow(data)),
                    ID1 = apply(data, 1, genotypeToNumeric, genotypeCol = 8),
                    ID2 = apply(data, 1, genotypeToNumeric, genotypeCol = 9))
  header <- c(
    "##fileformat=VCFv4.2",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
    "##FORMAT<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tID1\tID2"
  )

  writeLines(header, con = vcf_file)
  write.table(vcf, file = vcf_file, append = TRUE, quote = FALSE,
              sep ="\t", row.names = FALSE, col.names = FALSE)
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


