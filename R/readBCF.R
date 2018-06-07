#' Reads gentoype (0, 1, 2, NA coded) and allele depths from BCF. Optionally,
#' reads only from a specific region of the reference.
#'
#' @param path Path to VCF or BCF file
#' @param region A region (like 'Chr1:1-1000') to extract. `path` must be an
#'        indexed BCF if this option is used.
#' @param rowsAreSamples If TRUE, transpose GT and AD matricies so rows are
#'        samples and SNPs are columns. This is the opposite of a BCF.
#'
#' @export
readBCF_GTandAD = function(path, region=NULL, rowsAreSamples=T, minMAF=0.0, maxMissing=1) {
  if (is.null(region)) region = "";
  ret = readBCFQuery_(path, region)
  nsnp = length(ret$POS)
  # matrixify the matrices
  ret$GT = matrix(unlist(ret$GT, recursive = F), nrow=nsnp, byrow = T)
  ret$AD_ref = matrix(unlist(ret$AD_ref, recursive = F), nrow=nsnp, byrow = T)
  ret$AD_alt = matrix(unlist(ret$AD_alt, recursive = F), nrow=nsnp, byrow = T)
  # Adjust missing gentoypes: -1 is missing in C++ output
  ret$GT[ret$GT<0] = NA
  # Missing rate and AF
  N = ncol(ret$GT)
  ret$MissRate = rowSums(is.na(ret$GT))/N
  ret$AF = rowSums(ret$GT, na.rm = T)/(2*rowSums(!is.na(ret$GT)))
  # Recode 012 to be minor allele counts
  ret$GT_minor = ret$GT
  ret$GT_minor[ret$AF>0.5,] = 2 - ret$GT_minor[ret$AF>0.5,]
  ret$MAF = rowSums(ret$GT_minor, na.rm = T)/(2*rowSums(!is.na(ret$GT_minor)))
  ret$Ho = rowSums(ret$GT_minor==1, na.rm = T)/(2*rowSums(!is.na(ret$GT_minor)))
  # Postions are zero-based in BCF, we want 1-based
  ret$POS = ret$POS + 1
  snp.keep = ret$MAF >= minMAF & ret$MissRate <= maxMissing
  for (mat in c("GT", "AD_ref", "AD_alt", "GT_minor")) {
    # Keep good SNPs
    ret[[mat]] = ret[[mat]][snp.keep,]
    if (rowsAreSamples) {
      ret[[mat]] = t(ret[[mat]])
    }
  }
  for (vec in c("CHROM", "POS", "AF", "MissRate", "MAF", "Ho")) {
    ret[[vec]] = ret[[vec]][snp.keep]
  }
  ret
}

filter_GTandAD = function(GTAD, min.maf) {

}

#' Reads contig names and lengths from BCF file as a data frame.
#'
#' @param path Path to VCF or BCF file
#'
#' @export
readBCF_Contigs = function(path) {
  ret = readBCFContigs_(path)
  # Adjust missing gentoypes: -1 is missing in C++ output
  as.data.frame(ret)
}
