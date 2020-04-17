
# 抽平函数
phyRare = function(ps = ps,N = 3000){
  library(phyloseq)
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  otb = as.data.frame(t(vegan_otu(ps)))


  otb1 <- vegan::rrarefy(t(otb), N)

  ps = phyloseq(otu_table(as.matrix(otb1),taxa_are_rows = F),
                sample_data(ps),
                tax_table(ps)
                )
  ps
  return(ps)
}

