require(parallel)
#' Get kmers for a given k and a sequence
#' Does not consider before or after termini
slidingKmers <- function(seq, k){
  iter = 1:(nchar(seq)-k+1)
  kmers = substring(seq, iter, iter+k-1)
  
  kmers
}

#' Get kmers for a given k and a sequence
#' Does not consider before or after termini
slidingKmersInd <- function(start, end, k){
  len = end-start
  iter = start:(end-k+1)
  data.frame(start=iter, end=iter+k-1)
}


cacheKmers <- function(pfms, fasta, dir.out, cores=2){
  # Pre-generate negative data for different kmers
  all.kmers = sort(unique(sapply(pfms, ncol)))
  
  t = mclapply(X=all.kmers, mc.cores=cores, FUN=function(k){
    writeLines(sprintf('Generating kmers k=%s', k))
    # Get all k-mers
    neg = lapply(fasta, slidingKmers, k)
    neg = unlist(neg, F, F)
    # neg = table(neg)
    
    # Write out
    fout = file.path(dir.out, sprintf('%s.rds', k))
    writeLines(sprintf('Writing kmers k=%s: %s', k, fout))
    saveRDS(neg, fout)
  })
  
}

matchCutoffs <- function(pwms, kmer.dir, file.out){
  
  quant=c(0.75, 0.8, 0.85, 0.9, 0.95, 
          0.99, 0.995, 0.996, 0.997, 0.998, 0.999, 1.0)
  
  # Get all kmer vals we have
  k.vals = unique(as.character(sapply(pwms, ncol)))
  
  # Read in cached k-mer data
  writeLines('Reading in k-mer data from cache ...')
  kmer.files = file.path(kmer.dir, paste0(k.vals, '.rds'))
  kmer.data = lapply(kmer.files, readRDS)
  names(kmer.data) = gsub('.rds', '', basename(kmer.files))
 
  #mc.cores = cores, 
  cutoffs = lapply(X=1:length(pwms), FUN=function(i){
    
    pwm = pwms[[i]]
    # Get k value
    k = as.character(ncol(pwm))
    
    writeLines(sprintf('%s (%s/%s), k=%s', names(pwms)[i], i, length(pwms), k))
    # Load ORF k-mers if not loaded before
    kmer.freq = kmer.data[[k]]
    seqs = names(kmer.freq)
    
    # Compute scores (bottleneck!)
    writeLines('  scoring ...')
    sc = matchScore(seqs, pwm, F, both.strands = F)
    
    # Repeat duplicated things
    sc.mss = lapply(1:length(sc$mss), function(i) rep(sc$mss[[i]], kmer.freq[[i]]) )
    sc.mss = unlist(sc.mss, F, F)
    
    sc.css = lapply(1:length(sc$css), function(i) rep(sc$css[[i]], kmer.freq[[i]]) )
    sc.css = unlist(sc.css, F, F)
    
    li = list(mss=quantile(sc.mss, quant), css=quantile(sc.css, quant))
    
    print(li)
    
    li
  })
  names(cutoffs) = names(pwms)
  # cutoffs = do.call(rbind, cutoffs)
  
  saveRDS(cutoffs, file.out)
  return(cutoffs)
}

# setwd('~/Development/phd/tfbs_yeast/')
# source('jaspar.R')
# fa = readRDS('data/yeast_orf_trans_fasta.rds')
# pfms = read.jaspar.pfm('data/pfm_fungi.txt')[2]
# cacheKmers(pfms, fa, 'data/kmers_tmp', cores = 8)
