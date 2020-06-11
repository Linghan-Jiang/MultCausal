library(MatrixEQTL)
main <- function(cis=1e-6, tra=1e-6, dis=1e6){	
	cat('cis ', cis, '\n')
	cat('tra ', tra, '\n')
	cat('dis ', dis, '\n')
	useModel=modelLINEAR
	SNP_file_name = "data/mateqtl/arabi.snp"
	snps_location_file_name = "data/mateqtl/arabi.snploc"
	expression_file_name = "data/mateqtl/arabi.exp"
	gene_location_file_name = "data/mateqtl/arabi.geneloc"
	prefix = 'arabi'
	pvOutputThreshold_cis = cis
	pvOutputThreshold_tra = tra
	cisDist = dis
	output_file_name_cis = sprintf("data/mateqtl/%s.%s.%sDist.cis", prefix, pvOutputThreshold_cis, cisDist)
	output_file_name_tra = sprintf("data/mateqtl/%s.%s.%sDist.tra", prefix, pvOutputThreshold_tra, cisDist)
	
	## Load genotype data
	snps = SlicedData$new();
	snps$fileDelimiter = " ";      # the TAB character
	snps$fileOmitCharacters = "NA"; # denote missing values;
	snps$fileSkipRows = 1;          # one row of column labels
	snps$fileSkipColumns = 1;       # one column of row labels
	snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
	snps$LoadFile(SNP_file_name);

	## Load gene expression data
	gene = SlicedData$new();
	gene$fileDelimiter = " ";      # the TAB character
	gene$fileOmitCharacters = "NA"; # denote missing values;
	gene$fileSkipRows = 1;          # one row of column labels
	gene$fileSkipColumns = 1;       # one column of row labels
	gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
	gene$LoadFile(expression_file_name);

	## Run the analysis
	snpspos = read.table(snps_location_file_name, header = FALSE, stringsAsFactors = FALSE);
	genepos = read.table(gene_location_file_name, header = FALSE, stringsAsFactors = FALSE);

	me = Matrix_eQTL_main(
		snps = snps, 
		gene = gene, 
		output_file_name     = output_file_name_tra,
		pvOutputThreshold     = pvOutputThreshold_tra,
		useModel = useModel, 
		verbose = TRUE, 
		output_file_name.cis = output_file_name_cis,
		pvOutputThreshold.cis = pvOutputThreshold_cis,
		snpspos = snpspos, 
		genepos = genepos,
		cisDist = cisDist,
		pvalue.hist = FALSE,
		min.pv.by.genesnp = FALSE,
		noFDRsaveMemory = FALSE);

	## Results:
	cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
}

args = commandArgs(T)
tryCatch({
	cis = as.numeric(args[[1]])
	tra = as.numeric(args[[2]])
	dis = as.numeric(args[[3]])
	main(cis, tra, dis)
}, error = function(e){
	main()
})
