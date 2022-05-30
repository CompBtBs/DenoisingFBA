library(SAVER)

### set path
list_files=c('data/datasetGSE110949.csv','data/datasetE-GEOD-86618.csv','data/datasetGSE118056.csv')
list_genes=c('data/met_genes_datasetGSE110949.txt','data/met_genes_datasetE-GEOD-86618.txt','data/met_genes_datasetGSE118056.txt')
list_output=c('data/datasetGSE110949_saver_counts.csv','data/datasetE-GEOD-86618_saver_counts.csv','data/datasetGSE118056_saver_counts.csv')
i=0
for (i in 1:1){
  input_file=list_files[i]
  input_gene_file=list_genes[i]


  ###read data
  raw.data <- read.table(input_file, header = TRUE, 
                         skip = 0, row.names = 1,   
                         check.names = FALSE)
  cortex <- as.matrix(raw.data)
  ### to save computational time, the prediction is done only on the metabolic genes
  genes=read.table(input_gene_file,              # TXT data file indicated as string or full path to the file
             header = FALSE,    # Whether to display the header (TRUE) or not (FALSE)
             sep = "",          # Separator of the columns of the file
             dec = ".")   

  genes <- unlist(genes)
  genes.ind <- which(rownames(cortex) %in% genes)
  #######
  
  cortex.saver <- saver(cortex, 
                        pred.genes = genes.ind,
                        pred.genes.only = TRUE, 
                        estimates.only = TRUE,
                        ncores = 4,
                        size.factor = 1               #no normalization is required
                        )
  
  ######
  write.csv(cortex.saver,file=list_output[i])
}