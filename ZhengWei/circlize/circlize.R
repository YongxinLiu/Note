library(statnet)
library(circlize)

data = read.csv("SC.csv",header=T,row=1)
my.data=as.matrix(data)
rownames(my.data) =c("CCK", "CNPK", "GCCK", "GCNPK")

colnames(my.data) =c("Alphaproteobacteria","Betaproteobacteria","Gammaproteobacteria",
                      "Deltaproteobacteria","Acidobacteria","Actinobacteria",
                      "Bacteroidetes","Chloroflexi","Firmicutes",  
                      "Gemmatimonadetes","Planctomycetes","Thaumarchaeota" ,
                      "Verrucomicrobia","Ascomycota",  "Basidiomycota", 
                       "Zygomycota")


grid.col = NULL
grid.col[c("CCK", "CNPK", "GCCK", "GCNPK")] = c("chocolate", "chocolate", "chocolate", "chocolate")
grid.col[colnames(my.data)] = c("lavender", "khaki","mistyrose", 
                                "sienna1", "skyblue", "brown1", 
                                "gold", "maroon", "salmon", "moccasin",
                                "wheat","black","green","cyan","pink","orange")

circos.par(gap.degree = c(rep(2, nrow(my.data)-1), 10, rep(2, ncol(my.data)-1), 10),
           start.degree = 180)

chordDiagram(my.data,
             directional = TRUE,
             diffHeight = 0.06,
             grid.col = grid.col, 
             transparency = 0.5)

legend("right",pch=20,legend=colnames(my.data),
       col=grid.col[colnames(my.data)],bty="n",
       cex=1,pt.cex=3,border="black")
circos.clear()
