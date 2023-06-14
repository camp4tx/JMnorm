# Author: Guanjue Xiang
# Date: 2023-03-03
# Name: JMnorm.script.R
# Version: 1.0
# Description: This script is used to jointly normalize multiple epigenomic features between target and reference
# Usage: 
# source('JMnorm.script.R')
# 'ref_ave_sig' and 'target_sig' should be two input matrice (N-by-M matrix), where N is the number of cCREs and M is the number of chromatin features
# K = determineK(ref_ave_sig, sample_size=20000, add_sn=1)
# target_sig_JMnorm = JMnorm(target_sig, ref_ave_sig, K)

# Use Joint Multi-feature normalization (JMnorm) to do normalization across features matrix between Target and Reference
JMnorm = function(ct1_rep1_i_sig, file_ref_sig, K){
	library(dynamicTreeCut)
	# Initial round of PCA QTnorm to do normalization across features matrix between Target and Reference
	ct1_rep1_i_sig_qt = PCA_QTnorm_mat(ct1_rep1_i_sig, file_ref_sig)

	# Generate PCA signal matrix and get PCA rotation matrix based on reference signal matrix
	pca_ref0 = prcomp(file_ref_sig, center = F,scale. = F)
	# get PCA rotation matrix
	pca_ref0_rotation = pca_ref0$rotation

	# K-means clustering on reference signal matrix at PCA space
	# get the signal matrix at PCA space in reference
	file_ref_sig_pcs = pca_ref0$x

	# K-means clustering on reference signal matrix at PCA space, K equals to the number of clusters in reference hclust + dynamicTreeCut 
	ref_km = kmeans(file_ref_sig_pcs, K)

	# get the mean signal of each cluster
	ref_cluster_id = ref_km$cluster
	ref_sig_km_mean = c()
	for (i in 1:K){
		# the number of data point in the cluster is smaller than 2, then use the sample's signal vector as the cluster mean vector
		if (sum(ref_cluster_id==i)<2){
			ref_sig_km_mean = rbind(ref_sig_km_mean, file_ref_sig_pcs[ref_cluster_id==i,])
		} else {
			ref_sig_km_mean = rbind(ref_sig_km_mean, colMeans(file_ref_sig_pcs[ref_cluster_id==i,]))
		}
	}

	# Assign each sample to the closest cluster based on the distance between the target sample signal and the reference cluster mean
	# get the PCA QTnorm normalized	signal matrix at PCA space in target
	ct1_rep1_i_sig_pcs = ct1_rep1_i_sig_qt %*% pca_ref0_rotation

	# calculate the distance between each sample and the reference cluster mean at PCA space
	dist_to_ref_epistate_ct1_rep1 = pairwise_dist_mat((ct1_rep1_i_sig_pcs), (ref_sig_km_mean))
	# assign each sample to the closest cluster
	tar_cluster_id = apply(as.matrix(dist_to_ref_epistate_ct1_rep1), 1, function(x) which.min(as.numeric(x)))

	# Use PCA QTnorm to do normalization across features matrix between Target and Reference for each cluster
	# initialize the JMnorm normalized signal matrix at PCA space in target
	ct1_rep1_i_sig_JMnorm_pcs = ct1_rep1_i_sig_pcs
	# loop through each cluster
	for (ki in unique(tar_cluster_id)){
		# if the cluster has more than one sample
		if (sum(tar_cluster_id==ki)>1){
			# get the number of samples in the cluster in reference and target
			print(c(ki, sum(ref_cluster_id==ki), sum(tar_cluster_id==ki) ))
			# get the signal matrix at PCA space in reference and target in the cluster
			pca_ref0_ki_pcs = as.matrix(file_ref_sig_pcs[ref_cluster_id==ki, ])
			pca_tar1_ki_pcs = as.matrix(ct1_rep1_i_sig_pcs[tar_cluster_id==ki, ])

			# Use quantile normalization to do normalization across features matrix at PCA space between Target and Reference in the cluster 
			pca_tar1_ki_pcs_qt = pca_tar1_ki_pcs
			for (i in 1:dim(pca_tar1_ki_pcs)[2]){
				pca_tar1_ki_pcs_qt[,i] = QTnorm(pca_tar1_ki_pcs[,i], pca_ref0_ki_pcs[,i])
			}

			# update the normalized signal matrix at PCA space in target
			ct1_rep1_i_sig_JMnorm_pcs[tar_cluster_id==ki,] = pca_tar1_ki_pcs_qt
		}
	}

	# transform the normalized signal matrix at PCA space in target to the original space
	ct1_rep1_i_sig_JMnorm = ct1_rep1_i_sig_JMnorm_pcs %*% t(pca_ref0_rotation)
	return(ct1_rep1_i_sig_JMnorm)
}

# calculate the R2 value between two vectors
r2 = function(tar, ref){
	ssr1 = sum((tar-(ref))^2)
	ssr2 = sum((ref-(tar))^2)
	sst = sum((ref-mean(ref))^2)
	mse12 = 1 - ssr1/sst*0.5 - ssr2/sst*0.5
	return(mse12)
}

# quantile normalization for target_sig and ref_sig
QTnorm = function(target_sig,ref_sig){
	set.seed(2019)
	if (length(target_sig)==length(ref_sig)){
		# quantile normalization when target_sig is equal to ref_sig
		target_sig_QTnorm = sort(ref_sig)[rank(target_sig)]
	} else if (length(target_sig)>length(ref_sig)) {
		# quantile normalization when target_sig is larger than ref_sig
		ref_sig_s = sample(ref_sig, length(target_sig), replace=T)
		target_sig_QTnorm = sort(ref_sig_s)[rank(target_sig)]
	} else{
		# quantile normalization when target_sig is smaller than ref_sig
		ref_sig_s = sample(ref_sig, length(target_sig), replace=F)
		target_sig_QTnorm = sort(ref_sig_s)[rank(target_sig)]		
	}
	#target_sig_QTnorm = sort(ref_sig)[rank(target_sig)]
	target_sig_QTnorm[target_sig=0] = 0
	rm(target_sig)
	rm(ref_sig)
	return(target_sig_QTnorm)
}

# Use s3norm to normalize Signal to Noise Ratio (SNR) across features
normalize_snr = function(file_ref_sig0){
	# Use s3norm to normalize Signal to Noise Ratio (SNR) across features
	# Input: file_ref_sig0, a matrix of reference signal
	# Output: file_ref_sig, a matrix of normalized reference signal
	# get common pk based on quantiles
	file_ref_sig_topmean = apply(file_ref_sig0, 2, function(x) mean(x[x>=quantile(x, 0.99)]))
	file_ref_sig_colmean = apply(file_ref_sig0, 2, mean)
	# get common pk reference mean to normalize against
	ref_top_mean = median(file_ref_sig_topmean)
	# get common bg reference mean to normalize against
	ref_colmean = median(file_ref_sig_colmean)
	# get the Beta and Alpha in S3norm
	beta = (ref_top_mean-ref_colmean) / (file_ref_sig_topmean - file_ref_sig_colmean)
	alpha = ref_top_mean - beta * file_ref_sig_topmean

	# Use s3norm to normalize each feature against reference
	file_ref_sig = t(apply(file_ref_sig0, 1, function(x)  x * beta + alpha  ))
	# normalize each feature to have mean equals 1
	file_ref_sig = apply(file_ref_sig, 2, function(x) x-mean(x)+1)
	file_ref_sig[file_ref_sig<0] = 0
	return(file_ref_sig)
}

# Determine the K value for JMnorm based on reference signal matrix
determineK = function(refe_ave_sigmat, sample_size=20000, add_sn=1){
	library(dynamicTreeCut)
	# normalize signal to noise ratio (SNR) across features
	file_ref_sig = normalize_snr(refe_ave_sigmat)
	file_ref_sig[file_ref_sig<0] = 0

	# Cluster the reference signal matrix based on signal at PCA space
	set.seed(2019)
	# log2 transformation
	file_ref_sig_log2 = log2(file_ref_sig+add_sn)

	# if the number of rows is larger than 20000, randomly select 20000 rows to do PCA
	if (dim(file_ref_sig_log2)[1]>sample_size){
		file_ref_sig_log2_s = file_ref_sig_log2[sample(dim(file_ref_sig_log2)[1], sample_size),]
	}else{
		file_ref_sig_log2_s = file_ref_sig_log2
	}

	# PCA transformation
	pca_ref0 = prcomp(file_ref_sig_log2_s, center = F,scale. = F)

	# Get the distance matrix based on PCA space
	cCRE_dist_ref = dist(pca_ref0$x)

	# Cluster the reference signal matrix based on signal at PCA space by using Hclust and cutreeDynamic
	ref_km = hclust(cCRE_dist_ref, method='complete')
	ref_km_cluster = cutreeDynamic(ref_km, minClusterSize=round(dim(file_ref_sig_log2_s)[1]*0.01), distM=as.matrix(cCRE_dist_ref), method='hybrid', deepSplit=2)

	# print the number of data points in each cluster
	print('cutreeDynamic: ')
	print(table(ref_km_cluster))

	# Get the number of clusters as K value
	K = max(ref_km_cluster)
	return(K)
}

# Use PCA QTnorm to do normalization across features matrix between Target and Reference
PCA_QTnorm_mat = function(ct1_rep1_i_sig, file_ref_sig){
	# Generate PCA signal matrix and get PCA rotation matrix based on reference signal matrix
	pca_ref0 = prcomp(file_ref_sig, center = F,scale. = F)
	# get PCA rotation matrix
	pca_ref0_rotation = pca_ref0$rotation
	# get the signal matrix at PCA space in reference and target
	pca_ref0_pcs = pca_ref0$x
	pca_tar1_pcs = as.matrix(ct1_rep1_i_sig) %*% pca_ref0_rotation
	# quantile normalization at PCA space
	pca_tar1_pcs_qt = pca_tar1_pcs
	for (i in 1:dim(pca_tar1_pcs)[2]){
		pca_tar1_pcs_qt[,i] = QTnorm(pca_tar1_pcs[,i], pca_ref0_pcs[,i])
	}
	# transform back to original space
	tar1PCA_QT = pca_tar1_pcs_qt %*% t(pca_ref0_rotation)
	return(tar1PCA_QT)
}


# Compute pairwise Euclidean distance between two matrices
pairwise_dist_mat <- function(mat1, mat2) {
	# Calculate the row sums of squared elements for each matrix
	mat1_sumsq <- rowSums(mat1^2)
	mat2_sumsq <- rowSums(mat2^2)
	
	# Calculate the cross-product of mat1 and the transpose of mat2
	cross_prod <- mat1 %*% t(mat2)
	
	# Compute the pairwise Euclidean distance matrix using the outer product
	# of the row sums of squared elements and the cross-product
	# dist_mat[i, j] = sqrt(mat1_sumsq[i] + mat2_sumsq[j] - 2 * cross_prod[i, j])
	dist_sq = outer(mat1_sumsq, mat2_sumsq, "+") - 2 * cross_prod
	dist_sq[dist_sq<0]=0
	dist_mat <- sqrt(dist_sq)
	return(dist_mat)
}


