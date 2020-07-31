log <- file(snakemake@log[[1]], open="wt")

cdf_t = snakemake@params[["ct"]]
min_rep = snakemake@params[["mr"]]
min.p.mean = snakemake@params[["mm"]]
path_run_metatda = snakemake@params[["pm"]]
path_delta = snakemake@params[["path_delta"]]
path_out = snakemake@params[["path_out"]]



library(data.table)
library(distributions3)





get_rep_table <- function( file_path, rep){


  
  comparison <- data.table()
  
  
  for ( i in seq(1:rep)){
    
    print(i)
    
    path <- paste0(file_path, i, ".diff")
    file <- fread(path)
    file[, Rep:=i]
    comparison <- rbind(comparison, file)
  
  }
  
  
  colnames(comparison) <- c( "Gene", "Node", "Coord", "Strand", "Type", "Psi_A", "Psi_B", "DeltaPsi", "Probability", "Complexity", "Entropy", "V1", "Rep")
  comparison[ , V1:=NULL]
  
  comparison

}






cdf.beta <- function(mu, var, p) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  
  return(cdf( Beta(alpha, beta), p))
}







get_diff_nodes <- function (path, comp_name, reps, beta_t, min.p.mean, min_number_reps){
  
  
    print(comp_name)

    comp <- get_rep_table( paste0(path, comp_name, "_rep_") , reps)
    
    
    
    comp.stats <- comp[ , .(Psi_A.mean =mean(Psi_A) , Psi_B.mean =mean(Psi_B) ,
                                                                                        DeltaPsi.mean = mean(DeltaPsi), DeltaPsi.sd = sd(DeltaPsi),
                                                                                        Probability.mean=mean(Probability, na.rm=T), Probability.sd=sd(Probability, na.rm=T),
                                                                                        Probability.var = var(Probability, na.rm=T),  Number=.N),
                                                                                    by=c( "Gene", "Node", "Coord", "Strand", "Type")]
    
    
    comp.stats[ ,  cdf.beta:=cdf.beta( Probability.mean, Probability.var  , beta_t) ]
    
    
    comp.stats[ , diff:=(abs(DeltaPsi.mean)>=0.2 & cdf.beta  < 0.05 & Probability.mean>=min.p.mean & ! Type %in% c("TE", "TS") & Number > min_number_reps)  ]
    


}


snakepool_BetaDist <-function(beta_t, min.p.mean, min_number_reps, path_metadata, path_delta, out_dir){
  
  


run_metadata <-  fread(path_metadata)



for (i in 1:nrow(run_metadata)) {
  
  #print(run_metadata[i, Compare_ID])
  #print(run_metadata[i, Repeat])
  
  out <- get_diff_nodes(path_delta, run_metadata[i, Compare_ID], run_metadata[i, Repeat], beta_t, min.p.mean, min_number_reps )
  
  fwrite(out, file= paste0(out_dir, run_metadata[i, Compare_ID], ".txt") , append = FALSE, quote = "auto", sep = "\t",  row.names = FALSE, col.names = TRUE)

}

snakepool_BetaDist(cdf_t, min.p.mean, min_rep,  
                   path_run_metatda,
                   path_delta,
                   path_out)



