#' batch_sil
#'
#' Determine batch/bio effect using the silhouette coefficient (adopted from scone):
#' @param pca.data a list as created by 'prcomp', batch_sil.R needs $x: the principal components (PCs, correctly: the rotated data)
#' @param batch vector with the batch covariate (for each cell)
#' @param nPCs the number of principal components to use (default: 3)
#' @return The average silhouette width for all clusters. For batch effect, the smaller the better. For biological effect, the larger the better.
#' @export
batch_sil <- function(pca.data, batch, nPCs=3){
  # in scone, they use svd to compute principal components.
  # For now, we'll keep the PCA object created in prcomp to be consistent
  #proj <- svd(scale(t(expr),center = TRUE,scale = TRUE), nu = 3, nv =0)$u

  dd <- as.matrix(dist(pca.data$x[,1:nPCs]))
  score_sil <- summary(cluster::silhouette(as.numeric(batch), dd))$avg.width

  return(score_sil)
}


