
#---------纯R语言绘制稀释曲线
# #---候选指标
# all = c("observed" , "chao1"  , "diversity_inverse_simpson" , "diversity_gini_simpson",
#         "diversity_shannon"   ,   "diversity_fisher"   ,  "diversity_coverage"     ,    "evenness_camargo",
#         "evenness_pielou"    ,   "evenness_simpson"       ,    "evenness_evar" ,   "evenness_bulla",
#         "dominance_dbp"      ,  "dominance_dmn"        ,      "dominance_absolute"   ,      "dominance_relative",
#         "dominance_simpson"      ,    "dominance_core_abundance" ,  "dominance_gini"  ,           "rarity_log_modulo_skewness",
#         "rarity_low_abundance"   ,    "rarity_noncore_abundance",  "rarity_rare_abundance")



# 测试代码
# otu = otu
# map = map


# otu = read.delim("./otutab.txt",row.names = 1)
# map = read.delim("./metadata.tsv",row.names = 1)
# result1 = rareAlpha(otu = otu,map = map,group = "Group", start = 100,step = 100,method = "chao1")
# # 提取图形
# p =  result1[[1]]
# p
# # 提取数据
# data = result1[[2]]
# head(data)
# # 提取分组稀释曲线
# p =  result1[[3]]
# p
# # 提取分组稀释曲线
# p =  result1[[4]]
# p

# ps = readRDS("./ps_liu.rds")
# result = rareAlpha(ps = ps, start = 100,step = 1000,method = "observed")
# # 提取图形
# p =  result[[1]]
# p
# # 提取数据
# data = result[[2]]
# head(data)
# # 提取分组稀释曲线
# p =  result[[3]]
# p
# # 提取分组稀释曲线
# p =  result[[4]]
# p

# 稀释曲线函数
rareAlpha =function(otu = NULL,map = NULL,ps = NULL,method = "Richness",group = "Group", start = 100,step = 100){
  # 所需R包
  library(vegan)
  library(microbiome)
  library(tidyverse)
  # 搜所需函数
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
                  sample_data(ps)
    )
    ps
    return(ps)
  }
  
  # 数据导入
  if (is.null(otu)&is.null(map)) {
    ps = ps
    
  }else {
    #导入otu表格
    otu = as.matrix(otu)
    str(otu)
    # #导入注释文件
    # tax = as.matrix(tax)
    # taxa_names(tax)
    
    #导入分组文件
    colnames(map) = gsub(group,"AA", colnames(map))
    
    map$Group = map$AA
    map$Group = as.factor(map$Group )
    map$Group 
    # #导入进化树
    # tree = read.tree("./otus.tree")
    # tree
    ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
                   sample_data(map) 
                   # phy_tree(tree)
    )
  }
  
  
  # 全部指标
  all = c("observed" , "chao1"  , "diversity_inverse_simpson" , "diversity_gini_simpson",
          "diversity_shannon"   ,   "diversity_fisher"   ,  "diversity_coverage"     ,    "evenness_camargo",
          "evenness_pielou"    ,   "evenness_simpson"       ,    "evenness_evar" ,   "evenness_bulla",
          "dominance_dbp"      ,  "dominance_dmn"        ,      "dominance_absolute"   ,      "dominance_relative",
          "dominance_simpson"      ,    "dominance_core_abundance" ,  "dominance_gini"  ,           "rarity_log_modulo_skewness",
          "rarity_low_abundance"   ,    "rarity_noncore_abundance",  "rarity_rare_abundance")
  # 运行计算
  for (i in seq(start,max(sample_sums(ps)), by = step) ) {
    
    psRe = phyRare(ps = ps,N = i)
    
    vegan_otu <-  function(physeq){
      OTU <-  otu_table(physeq)
      if(taxa_are_rows(OTU)){
        OTU <-  t(OTU)
      }
      return(as(OTU,"matrix"))
    }
    
    if (method == "Richness") {
      count = as.data.frame(t(vegan_otu(psRe)))
      head(count)
      
      # # 加载VEGAN包
      # alpha=diversity(count, "shannon")
      
      x = t(count) ##转置，行为样本，列为OTU
      head(x)
      est <- estimateR(x)
      est <- estimateR(x)
      index <- est[1, ]
      
    }
    
    if (method %in% c("ACE")) {
      ap_phy = estimate_richness(psRe, measures =method)
      head(ap_phy)
      index = ap_phy$ACE
      
    }
    
    if (method %in% all) {
      alp_mic = alpha(psRe,index=method)
      head(alp_mic)
      index= alp_mic[,1]
      
    }
    
   
    
    
    
    tab = data.frame(ID = names(sample_sums(psRe)))
    #-得到多样性的列
    tab = cbind(tab,i,index)
    head(tab)
    if (i == start) {
      result = tab
    }
    
    if (i != start) {
      result = rbind(result,tab)
    }
    
    
  }
 
   # ---------稀释结果调整
  for (ii in 1:length(sample_sums(ps))) {
    result$i[result$i > sample_sums(ps)[ii][[1]]]
    df_filter<- filter(result, ID ==names(sample_sums(ps)[ii]) &i > sample_sums(ps)[ii][[1]])
    result$index
    result$index[result$i>sample_sums(ps)[ii][[1]]]
    a = result$i>sample_sums(ps)[ii][[1]]
    a[a == FALSE] = "a"
    b = result$ID == names(sample_sums(ps)[ii])
    b[b == FALSE] = "b"
    result$index[a== b] = NA
  }
  
  map = as.data.frame(sample_data(ps))
  result$Group = map$Group
  ## 绘制稀释曲线
  
  library(ggplot2)
  
  p = ggplot(data= result,aes(x = i,y = index,group = ID,colour = Group)) +
    # geom_point() +
    # geom_line() +
    geom_smooth(span = 0.7, se = FALSE, method = "loess") +
    labs(x= "",
         y=method,
         title="") +theme_bw()+
    
    #scale_y_continuous(expand = c(0,0))+
    # geom_hline(aes(yintercept=0), colour="black", linetype=2) +
    # geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
    # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
    theme(
      
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      
      plot.title = element_text(vjust = -8.5,hjust = 0.1),
      axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
      axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
      axis.text = element_text(size = 20,face = "bold"),
      axis.text.x = element_text(colour = "black",size = 14),
      axis.text.y = element_text(colour = "black",size = 14),
      legend.text = element_text(size = 15,face = "bold")
      #legend.position = "none"#是否删除图例
      
    ) 
  
  p
  
  # 数据再分析
  data = result
  # ---分组求均值
  groups<- group_by(data, Group,i)
  data2 <- summarise(groups , mean(index), sd(index))
  # head(data2)
  colnames(data2) = c(colnames(data2)[1:2],"mean","sd")
  # 再出图
  p2 = ggplot(data= data2,aes(x = i,y = mean,colour = Group)) +
    # geom_point() +
    # geom_line() +
    geom_smooth(span = 0.7,se = FALSE, method = "loess") +
    # geom_errorbar(data = data3,aes(ymin=mean - sd, ymax=mean + sd), colour="black", width=.1)+
    labs(x= "",
         y=method,
         title="") +theme_bw()+
    
    #scale_y_continuous(expand = c(0,0))+
    # geom_hline(aes(yintercept=0), colour="black", linetype=2) +
    # geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
    # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
    theme(
      
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      
      plot.title = element_text(vjust = -8.5,hjust = 0.1),
      axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
      axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
      axis.text = element_text(size = 20,face = "bold"),
      axis.text.x = element_text(colour = "black",size = 14),
      axis.text.y = element_text(colour = "black",size = 14),
      legend.text = element_text(size = 15,face = "bold")
      #legend.position = "none"#是否删除图例
      
    ) 
  
  p2
  # 再出图
  p4 = ggplot(data= data2,aes(x = i,y = mean,colour = Group)) +
    # geom_point() +
    # geom_line() +
    # geom_smooth(span = 0.7,se = FALSE, method = "loess") +
    geom_errorbar(data = data2,aes(ymin=mean - sd, ymax=mean + sd,colour = Group),alpha = 0.4, width=.1)+
    labs(x= "",
         y=method,
         title="") +theme_bw()+
    
    #scale_y_continuous(expand = c(0,0))+
    # geom_hline(aes(yintercept=0), colour="black", linetype=2) +
    # geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
    # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
    theme(
      
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      
      plot.title = element_text(vjust = -8.5,hjust = 0.1),
      axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
      axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
      axis.text = element_text(size = 20,face = "bold"),
      axis.text.x = element_text(colour = "black",size = 14),
      axis.text.y = element_text(colour = "black",size = 14),
      legend.text = element_text(size = 15,face = "bold")
      #legend.position = "none"#是否删除图例
      
    )
  
  p4
  return(list(p,table = result,p2,p4))
}



