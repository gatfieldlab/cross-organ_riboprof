#function to get dissimilarity matrix between the 4 samples (RPF liver, RPF kidney, RNA liver, RNA kidney).Only for 1 gene!
get_dist_matrix <- function(gene, norm_mean_data=rpkm.mean.log.comb, metric='euclidean'){
  norm_matrix <- t(data.frame(Liver_RNA = norm_mean_data[gene, 1:12], Liver_RPF = norm_mean_data[gene, 13:24], 
                              Kidney_RNA = norm_mean_data[gene, 25:36], Kidney_RPF = norm_mean_data[gene, 37:48]))
  dist(norm_matrix, method=metric)
}

#FUNCION To get just expression matrix
get_expression_matrix <- function(genes, norm_mean_data=rpkm.mean.comb){
  norm_matrix <- data.frame(row.names = c('Liver_RNA','Liver_RPF','Kidney_RNA','Kidney_RPF'))
  for(gene in genes){
    norm_matrix <- cbind(norm_matrix, t(data.frame(Liver_RNA = norm_mean_data[gene, 1:12], Liver_RPF = norm_mean_data[gene, 13:24], 
                                                   Kidney_RNA = norm_mean_data[gene, 25:36], Kidney_RPF = norm_mean_data[gene, 37:48])))
  }
  norm_matrix
}

#function to draw cluster trees for individual genes.
get_cluster_tree <- function(genes, norm_mean_data=rpkm.mean.comb, method='average', save=FALSE){
  for (gene in genes) {
    if (save){
      pdf(paste(gene, annot_gene_names[gene,], method, 'withAC_clustering.pdf', sep='_'))
      plot(agnes(get_dist_matrix(gene, norm_mean_data=norm_mean_data),diss=TRUE,method=method), xlab=paste(gene, annot_gene_names[gene,], sep='_'))
      #plot(hclust(get_dist_matrix(gene), method=method), xlab=paste(gene, annot_gene_names[gene,], sep='_'))
      dev.off()
    } else {
      #plot(hclust(get_dist_matrix(gene), method=method), xlab=paste(gene, annot_gene_names[gene,], sep='_'))
      plot(agnes(get_dist_matrix(gene, norm_mean_data=norm_mean_data),diss=TRUE,method=method), xlab=paste(gene, annot_gene_names[gene,], sep='_'))
    }}
}

get_cluster_tree(genes = core_clock2$gid)


#weighted average of dissimilarities matrices
library(analogue) #for fuse function
get_mean_dist <- function(genes, norm_data) {
  genesdist <- lapply(genes, function(x){get_dist_matrix(x,norm_mean_data=norm_data)})
  do.call(fuse, genesdist)
}

plot(hclust(get_mean_dist(genes = core_clock2$gid, norm_data=rpkm.mean.comb),method='average')) #for raw rpkm data

