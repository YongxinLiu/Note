# 制作树图格式数据 Format feature table and taxonomy into maptree data
#
# This is the function named 'format2maptree'
# which format feature table and taxonomy into maptree data, and reture a data object
#
#' @title Format feature table and taxonomy into maptree data
#' @description Input taxonomy composition, and feature table. Then select top N high abundance features.Finally, return a maptree(ggraph) object.
#' @param otu OTU/ASV table;
#' @param tax taxonomy annotation, include feature ID and 7 levels;
#' @param N Top N features to show, default 200;
#' @details
#' By default, returns top 8 taxonomy and group mean stackplot
#' The available style include the following:
#' \itemize{
#' \item{group: group mean circlize}
#' \item{sample: each sample circlize}
#' }
#' @return ggplot2 object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
#' Root microbiota shift in rice correlates with resident time in the field and developmental stage.
#' Sci China Life Sci 61, DOI: \url{https://doi.org/10.1007/s11427-018-9284-4}
#'
#' @seealso tax_maptree
#' @examples
#' # Input feature table, taxonomy and Top N features, and format into mapdata
#' mapdata = format2maptree(otutab, taxonomy, 200)
#' # Add mean abundance size and phylum color for maptree
#' mapadd = tax_maptree(mapdata)
#' # Saving and plotting maptree
#' (p = mapadd[[1]])
#' @export


format2maptree = function(otu = otutab, tax = taxonomy, N = 200){
  ##-------------------------------------------依赖包-----------
  suppressWarnings(suppressMessages(library(igraph)))
  suppressWarnings(suppressMessages(library(ggraph)))
  suppressWarnings(suppressMessages(library(tidyverse)))
  suppressWarnings(suppressMessages(library(viridis)))
  suppressWarnings(suppressMessages(library(data.tree)))
  suppressWarnings(suppressMessages(library(phyloseq)))
  suppressWarnings(suppressMessages(library(ggtree)))
  suppressWarnings(suppressMessages(library(RColorBrewer)))

  # 导入功能函数
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  vegan_tax <-  function(physeq){
    tax <-  tax_table(physeq)

    return(as(tax,"matrix"))
  }
  ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE),
                 tax_table(tax)
  )
  #相对丰富转换
  ps1_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela
  # 提取otu表格
  otu_table = as.data.frame(t(vegan_otu(ps1_rela)))
  otu_table$mean = rowMeans(otu_table)
  otu_table$ID = row.names(otu_table)
  head(otu_table)
  #按照从大到小排序
  otu_table<- arrange(otu_table, desc(mean))
  subtab = head(otu_table,N)
  head(subtab)
  row.names(subtab) =subtab$ID
  subtab = subtab[,1:(dim(subtab)[2]-2)]
  subtab = as.matrix(subtab)

  #对phyloseq取子集
  ps_sub <- phyloseq(otu_table(subtab, taxa_are_rows=TRUE),
                     tax_table(tax))
  ps_sub
  #--------------------整理tax文件-----------------
  #注意我们本次制作的是OTU水平的maptree
  tax = as.data.frame(vegan_tax(ps_sub))
  #添加水平识别标示
  tax$Kingdom = paste("D_",tax$Kingdom,sep = "--")
  tax$Phylum = paste(tax$Kingdom,"P_",tax$Phylum,sep = "--")
  tax$Class = paste(tax$Phylum,"C_",tax$Class,sep = "--")
  tax$Order = paste(tax$Class,"O_",tax$Order,sep = "--")
  tax$Family = paste(tax$Order,"F_",tax$Family ,sep = "--")
  tax$Genus = paste(tax$Family,"G_",tax$Genus,sep = "--")
  tax$Species = paste(tax$Genus,"S_",tax$Species,sep = "--")
  tax$OTU = paste("ASV",row.names(tax),sep = "--")
  # # 对不同分类水平的tax未知物种进行区分
  tax_table(ps_sub) = as.matrix(tax)
  head(tax)

  #--------------------构造边文件--------------------
  #提取分组信息
  Desep_group <- as.character(colnames(tax))
  #按照form--to格式来排布数据
  edge_t = tax[c(Desep_group[1:2])]
  dim(edge_t)
  edge_t <-  dplyr::distinct(edge_t)
  dim(edge_t)
  for (i in 2:7) {
    result = tax[c(Desep_group[i:(i+1)])]
    colnames(result) = colnames(edge_t)
    result <-  dplyr::distinct(result)
    edge_t = rbind(edge_t,result)
  }
  dim(edge_t)
  colnames(edge_t) = c("from","to")
  if (length(unique(tax$Kingdom)) >1) {
    edge_t$from = as.character(edge_t$from)
    edge_t$to = as.character(edge_t$to)
    buc = data.frame(from = c("king_up","king_up"),to =c("D_Archaea","D_Bacteria") )
    row.names(buc) = c("sp_K","sp_k1")
    deg = rbind(edge_t,buc)
  }
  deg = edge_t
  #---------------------------构造基本节点信息文件
  #构造tree用于提取分级结构
  tree <- FromDataFrameNetwork(deg)
  vertices_t  <-  data.frame(
    name = unique(c(as.character(deg$from), as.character(deg$to)))
  )
  # Then I can easily get the level of each node, and add it to the initial data frame:
  mylevels <- data.frame( name=tree$Get('name'), level=tree$Get("level") )
  vertices_t <- vertices_t %>%
    left_join(., mylevels, by=c("name"="name"))
  #这里提取名称信息
  AA = rep("A",length(vertices_t$name))
  for (i in 1:length(vertices_t$name)) {
    AA [i] = strsplit(basename(as.character(vertices_t$name)), "--")[[i]][length(strsplit(basename(as.character(vertices_t$name)), "--")[[i]])]
  }
  vertices_t$shortName<- AA
  # #设置level为2的定义为标签，其他偶读都省略掉，避免杂乱无章
  vertices_t <- vertices_t %>%
    mutate(new_label=ifelse(level==2, shortName, NA))
  row.names(vertices_t) = vertices_t$name

  #----------------------------------构造网络文件------------------------------
  mygraph <- graph_from_data_frame( deg, vertices= vertices_t )

  p = ggraph(mygraph, layout = 'circlepack', sort.by = NULL, direction = "out") +
    geom_node_circle(aes(fill = as.factor(depth),color = as.factor(depth) ) ) +
    scale_fill_manual(values= c("white","white","white","white","white","white","white","yellow") ) +
    scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black", "5"="black", "6"="black", "7"="black") ) +
    geom_node_text( aes(label=new_label), size=6,repel = TRUE) +
    # geom_node_label( aes(label=new_label), size=3,repel = TRUE) +
    theme_void() +
    theme( legend.position="FALSE",plot.margin = unit(rep(0,4), "cm"))#
  # p

  return(list(p,deg,vertices_t,ps_sub,mygraph))
}

