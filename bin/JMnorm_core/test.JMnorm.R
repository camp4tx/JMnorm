# Author: Guanjue Xiang
# Date: 2023-03-03
# Description: This script is used to jointly normalize multiple epigenomic features between target and reference

library(readr)
library(pheatmap)

# Source the JMnorm script 
source('/mnt/ebs/VISION_hg38_wg/JMnorm.script.R')
# Set the working directory
setwd('/mnt/epimap/cCRE_signal_mat/VISION')


# get cell-type list
ct_list = c('TCD8','TCD4','NK','B','ERY','MK','MONp','MONc','NEU')
# get feature list
feature_list = c("ATAC", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
# added 1 to avoid 0 in log transformation
add_sn = 1


# read refe_ave_sigmat
file_ref = as.data.frame(read_table('ref.raw_sigmat.rep1_2.txt', col_names = F))
file_ref_sig0 = file_ref[,-1]
# log2 transformation
file_ref_sig0_log2 = log2(file_ref_sig0+add_sn)


# Determine the K value for JMnorm based on reference signal matrix
K = determineK(file_ref_sig0, sample_size=20000, add_sn=add_sn)


# JMnorm for each cell-type ct_i
ct_i = ct_list[1]
# read target signal matrix
target_ct_rep1_i = as.data.frame(read_table(paste0(ct_i,'.raw_sigmat.rep1.txt'), col_names = F))
# add noise to avoid many 0s
target_ct_rep1_i_sig0 = target_ct_rep1_i[,-1]+matrix(runif(prod(dim(target_ct_rep1_i)), min = -0.00001, max = 0.0001), nrow = dim(target_ct_rep1_i)[1])
# avoid negative values
target_ct_rep1_i_sig0 = target_ct_rep1_i_sig0-min(target_ct_rep1_i_sig0)
# normalize signal to noise ratio (SNR) across features
target_ct_rep1_i_sig = normalize_snr(target_ct_rep1_i_sig0)
# log2 transformation
target_ct_rep1_i_sig_log2 = log2(target_ct_rep1_i_sig+add_sn)
# JMnorm
target_ct_rep1_i_sig_JMnorm = JMnorm(target_ct_rep1_i_sig_log2, file_ref_sig0_log2, K)

# add colnames
colnames(file_ref_sig0_log2) = feature_list
# correlation between features in reference cell-type
file_ref_sig0_log2_cor = cor(file_ref_sig0_log2)
# add colnames
colnames(target_ct_rep1_i_sig_JMnorm) = feature_list
# correlation between features in target cell-type before JMnorm
target_ct_rep1_i_sig_cor = cor(target_ct_rep1_i_sig_log2)
# correlation between features in target cell-type after JMnorm
target_ct_rep1_i_sig_JMnorm_cor = cor(target_ct_rep1_i_sig_JMnorm)

# calculate the R2 between the correlation matrix of target cell-type and reference cell-type
r2_after = r2(as.numeric(target_ct_rep1_i_sig_JMnorm_cor), as.numeric(file_ref_sig0_log2_cor))
r2_before = r2(as.numeric(target_ct_rep1_i_sig_cor), as.numeric(file_ref_sig0_log2_cor))
print(paste0('R2 of correlation matrix between target reference before JMnorm: ', r2_before))
print(paste0('R2 of correlation matrix between target reference after JMnorm: ', r2_after))

# plot the correlation matrix of target cell-type after JMnorm
breaksList = seq(-1,1, by = 0.001)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
pheatmap(target_ct_rep1_i_sig_JMnorm_cor, color=my_colorbar, breaks = breaksList, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T)
# plot the correlation matrix of target cell-type before JMnorm
pheatmap(target_ct_rep1_i_sig_cor, color=my_colorbar, breaks = breaksList, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T)
# plot the correlation matrix of reference cell-type
pheatmap(file_ref_sig0_log2_cor, color=my_colorbar, breaks = breaksList, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T)
