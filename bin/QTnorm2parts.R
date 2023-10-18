
# get z-score
get_z = function(x){
	x_notop = x[x<=quantile(x, 0.99)]
        xz = (x - mean(x_notop)) / sd(x_notop)
        return(xz)
}

# get fdr adjusted p-value
get_fdr = function(x){
        z = get_z(x)
        zp = pnorm(-abs(z))
        zpfdr = p.adjust(zp)
        return(zpfdr)
}

QTnorm = function(target_sig,ref_sig){
	set.seed(2019)
	if (length(target_sig)==length(ref_sig)){
		target_sig_qtnorm = sort(ref_sig)[rank(target_sig)]
	} else if (length(target_sig)>length(ref_sig)) {
		ref_sig_s = sample(ref_sig, length(target_sig), replace=T)
		target_sig_qtnorm = sort(ref_sig_s)[rank(target_sig)]
	} else{
		ref_sig_s = sample(ref_sig, length(target_sig), replace=F)
		target_sig_qtnorm = sort(ref_sig_s)[rank(target_sig)]		
	}
  target_sig_qtnorm[target_sig=0] = 0
  rm(target_sig)
  rm(ref_sig)
  return(target_sig_qtnorm)
}

QTnorm_2part = function(tar_sig, ref_sig, thresh){
	# sampling to match the length of two vectors
	if (length(tar_sig)>length(ref_sig)){
		ref_sig_s = sample(ref_sig, length(tar_sig), replace=T)
	} else {
		ref_sig_s = sample(ref_sig, length(tar_sig), replace=F)
	}
	# get binary vector for peak and background
	tar_sig_peak_binary = get_fdr(tar_sig)<thresh
	ref_sig_peak_binary = get_fdr(ref_sig)<thresh
	# peak
	tar_sig_pk = tar_sig[tar_sig_peak_binary]
	ref_sig_pk = ref_sig[ref_sig_peak_binary]
	# background
	tar_sig_bg = tar_sig[!tar_sig_peak_binary]
	ref_sig_bg = ref_sig[!ref_sig_peak_binary]
	# QTnorm peak & background
	tar_sig_pk_qt = QTnorm(tar_sig_pk, ref_sig_pk)
	tar_sig_bg_qt = QTnorm(tar_sig_bg, ref_sig_bg)
	# combine peak & background
	tar_sig_qt_2part = tar_sig
	tar_sig_qt_2part[tar_sig_peak_binary] = tar_sig_pk_qt
	tar_sig_qt_2part[!tar_sig_peak_binary] = tar_sig_bg_qt
	return(tar_sig_qt_2part)
}


