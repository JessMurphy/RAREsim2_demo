
# function to convert the haplotypes into a genotype matrix
make_geno = function(hap){
  
  # create an empty genotype matrix
  geno = matrix(0, nrow(hap), ncol(hap)/2)
  
  # sum up the number of alleles in adjacent haplotypes (2 haplotypes per person)
  for (j in 1:(ncol(hap)/2)){
    geno[,j] = hap[,2*j] + hap[,2*j-1]
  }
  geno = as.data.frame(geno)
  
  return(geno)
}

# function to calculate the allele counts/frequencies
calc_allele_freqs = function(geno, n, sum=NULL, source=NULL) {
  
  #counts = data.frame(count = rowSums(geno)) %>%
  # mutate(mac = ifelse(count>n, 2*n-count, count)) %>%
  #  mutate(maf = mac/(2*n))

  if (!is.null(sum)) {
    counts = data.frame(ac = sum) %>% mutate(af = ac/(2*n))
  } else {
    counts = data.frame(ac = rowSums(geno)) %>% mutate(af = ac/(2*n))
  }
  
  if (!is.null(source)) {counts$source = source}
  
  return(counts)
}

# function to create a dataframe with a line for each variant instead of just counts 
# (necessary for ProxECAT v2)
make_long = function(counts, leg, case, group){
  
  # add information to the counts
  temp = counts %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case=case, group=group)
  
  # remove the monomorphic variants 
  temp2 = temp %>% dplyr::filter(ac!=0)
  
  # repeat each variant mac times
  out = data.frame(lapply(temp2, rep, temp2$ac)) %>% select(-ac, -af)
  
  return(out)
}

# function to merge the case datasets for power and type I error calculations
merge_cases = function(cases.power, cases.t1e, leg, genes.power) {
  
  # add row number and gene column to each hap
  hap.power = cases.power %>% mutate(row = leg$row, gene = leg$gene)
  hap.t1e = cases.t1e %>% mutate(row = leg$row, gene = leg$gene)
  
  # subset haps to the necessary genes
  power.gene = subset(hap.power, gene %in% genes.power) 
  t1e.gene = subset(hap.t1e, !(gene %in% genes.power))
  
  # make sure the column names are the same
  names(t1e.gene) = names(power.gene)
  
  # merge the two case haps
  hap.out = rbind(power.gene, t1e.gene)
  
  # order the merged hap file by row number
  hap.out = hap.out[order(hap.out$row),]
  
  # remove the row number and gene columns
  hap.out = subset(hap.out, select = -c(row, gene))
  
  return(hap.out)
}

# function to flip values for a specified file at variants with an AF >= 1-maf 
# (used in the flip_data function below)
flip_file = function(file_to_flip, flip, file_type, N=NULL) {
  
  # Make copy of original file
  file2 = file_to_flip
  
  if (file_type == "leg") {
    
    # Flip ref allele in file2 with alt allele in file
    file2[flip, "a0"] <- file_to_flip[flip, "a1"]
    
    # Flip alt allele in leg2 with ref allele in leg
    file2[flip, "a1"] <- file_to_flip[flip, "a0"]
    
    return(file2)
    
  } else if (file_type == "geno") {
    
    # Flip the alternate allele counts at the relevant variants
    file2[flip, ] <- 2-file_to_flip[flip, ]
    
    return(file2)
    
  } else if (file_type == "count") {
    
    # Flip the ACs at the variants that need to be flipped
    file2[flip, "ac"] <- (2*N)-file_to_flip[flip, "ac"]
    
    # Flip the AFs at the variants that need to be flipped
    file2[flip, "af"] <- 1-file_to_flip[flip, "af"]
    
    return(file2)
    
  } else {
    stop("ERROR: 'file_type' must be a string of either 'leg', 'geno', or 'count'")
  }
  
}


# function to flip the relevant datasets at variants where AF >= 1-maf for each possible scenario
flip_data = function(leg, flip, geno_case, count_case, Ncase, cntrl, geno_ic=NULL, count_ic=NULL, Nic=NULL, geno_cc=NULL, count_cc=NULL, Ncc=NULL) {
    
    # Create new leg file
    leg2 = flip_file(leg, flip, file_type="leg")
    
    # Create new case data files
    geno_case2 = flip_file(geno_case, flip, file_type="geno")
    count_case2 = flip_file(count_case, flip, file_type="count", N=Ncase)
    
    if (cntrl == "int") {
      
      # Update geno files
      geno_ic2 = flip_file(geno_ic, flip, file_type="geno")
      
      # Recalculate ac/af 
      count_ic2 = flip_file(count_ic, flip, file_type="count", N=Nic)
      
      # Return all changed files, note some may be NULL
      return(list(leg=leg2, geno.case=geno_case2, geno.ic=geno_ic2, count.case=count_case2, count.ic=count_ic2))
      
    } else if (cntrl == "ext") {
      
      # Update geno files
      geno_cc2 = flip_file(geno_cc, flip, file_type="geno")
      
      # Recalculate ac/af 
      count_cc2 = flip_file(count_cc, flip, file_type="count", N=Ncc)
      
      # Return all changed files, note some may be NULL
      return(list(leg=leg2, geno.case=geno_case2, geno.cc=geno_cc2, count.case=count_case2, count.cc=count_cc2))
      
    } else if (cntrl == "all") {
      
      # Update geno files
      geno_ic2 = flip_file(geno_ic, flip, file_type="geno")
      geno_cc2 = flip_file(geno_cc, flip, file_type="geno")
      
      # Recalculate ac/af 
      count_ic2 = flip_file(count_ic, flip, file_type="count")
      count_cc2 = flip_file(count_cc, flip, file_type="count")
      
      # Return all changed files, note some may be NULL
      return(list(leg=leg2, geno.case=geno_case2, geno.ic=geno_ic2, geno.cc=geno_cc2, count.case=count_case2, count.ic=count_ic2, count.cc=count_cc2))
    }
}

# function for formatting data and running statistical test for ProxECAT by gene
prox_gene_data_prep = function(data_prox, common) {
  
  # count the number of fun and syn alleles by case status
  # (need .drop param so it still creates a group even if AC is 0)
  counts_gene = data_prox %>% count(gene, case, fun, .drop = FALSE)
  
  # convert the counts to wide format
  counts_wide = tidyr::pivot_wider(counts_gene, names_from=c(case, fun), values_from=n,
                                   values_fill=0, names_sep="_") %>% 
    mutate(case_ratio = case_fun/case_syn, control_ratio = control_fun/control_syn)
  
  # calculate medians (set na.rm to TRUE to avoid median ratio being NA if there's a divide by 0 problem)
  median_case_ratio = median(counts_wide$case_ratio, na.rm = TRUE)
  median_control_ratio = median(counts_wide$control_ratio, na.rm = TRUE)
  
  # calculate the necessary variables for proxecat-weighted and 
  # call proxecat if there are enough functional/synonymous variants
  counts_wide2 = counts_wide %>% mutate(case_fun_w = case_fun/median_case_ratio,
                                        control_fun_w = control_fun/median_control_ratio) %>%
  mutate(prox = ifelse(case_fun + control_fun < 5 | case_syn + control_syn < 5, NA,
                       proxecat(case_fun, case_syn, control_fun, control_syn)$p.value),
         prox_w = ifelse(case_fun_w + control_fun_w < 5 | case_syn + control_syn < 5, NA,
                         proxecat(case_fun_w, case_syn, control_fun_w, control_syn)$p.value))
  
  return(counts_wide2)
}

# function for formatting data and running statistical test for LogProx
logprox_gene_data_prep = function(data_logprox, current_gene, all=F, multiple=F, adjusted=F) {
  
  # filter data by gene
  data_gene = data_logprox %>% filter(gene==current_gene)
  
  # count the number of fun and syn alleles by case status
  # need .drop param so it still creates a group even if AC is 0
  counts_data_gene = data_gene %>% count(case, fun, .drop = FALSE)
  
  # fit the logprox model
  if (multiple){
    
    prox = tryCatch(ifelse((counts_data_gene$n[1] + counts_data_gene$n[3] < 5) | (counts_data_gene$n[2] + counts_data_gene$n[4] < 5),
                           NA, summary(glm(fun ~ case + source, data=data_gene, family="binomial"))$coefficients[2,4]), error = function(e) NA)
        
  } else if (all & adjusted){
    
    prox = tryCatch(ifelse((counts_data_gene$n[1] + counts_data_gene$n[3] < 5) | (counts_data_gene$n[2] + counts_data_gene$n[4] < 5),
                           NA, summary(glm(fun ~ case + group, data=data_gene, family="binomial"))$coefficients[2,4]), error = function(e) NA)
  
  } else{
    
    prox = tryCatch(ifelse((counts_data_gene$n[1] + counts_data_gene$n[3] < 5) | (counts_data_gene$n[2] + counts_data_gene$n[4] < 5),
                           NA, summary(glm(fun ~ case, data=data_gene, family="binomial"))$coefficients[2,4]), error = function(e) NA)
  }
    
  return(prox) 
}

# define function for calculating power
my.power = function(values, alpha=0.05){
  values2 = values[!is.na(values)]
  sig = which(as.numeric(values2) <= alpha)
  out = length(sig)/length(values2)
  return(out)
}
                                