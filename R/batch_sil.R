#steps to determine batch/bio effect in scone:
#' @export
batch_sil <- function(pca.data, batch, nPCs=3){
  # in scone, they use svd to compute principal components.
  # For now, we'll keep the PCA object created in prcomp to be consistent
  #proj <- svd(scale(t(expr),center = TRUE,scale = TRUE), nu = 3, nv =0)$u

  dd <- as.matrix(dist(pca.data$x[,1:nPCs]))
  score_sil <- summary(cluster::silhouette(as.numeric(batch), dd))$avg.width

  return(score_sil)
}


