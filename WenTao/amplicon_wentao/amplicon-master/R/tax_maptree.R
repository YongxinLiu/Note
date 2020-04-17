# 物种组成树图 Maptree of taxonomy composition
#
# This is the function named 'tax_maptree'
# which draw maptree, and reture a ggplot2 object
#
#' @title Plotting maptree of taxonomy
#' @description Visualize mapdata. Finally, return a ggplot2 object.
#' @param mapdata format2maptree result
#' @details
#' By default, returns  mapadd[[1]] is ggplot object
#' The available style include the following:
#' \itemize{
#' \item{phylum: color by phylum}
#' \item{diff: color by diff}
#' }
#' @return ggplot2 object.
#' @author Contact: Tao Wen, Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Jingying Zhang, Yong-Xin Liu, et.al. (2019).
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology 37, 676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @seealso data_to_maptree
#' @examples
#' # Input feature table, taxonomy and Top N features, and format into mapdata
#' mapdata = format2maptree(otutab, taxonomy, 200)
#' # Add mean abundance size and phylum color for maptree
#' mapadd = tax_maptree(mapdata)
#' # Saving and plotting maptree
#' (p = mapadd[[1]])
#' @export

##----注释和可视化树图Annotate and plot maptree----
tax_maptree = function(x){
  ps_sub = mapdata[[4]]
  vertices_t = mapdata[[3]]
  deg = mapdata[[2]]
  head(vertices_t)

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

  tax_table = as.data.frame(vegan_tax(ps_sub))
  otu_table = as.data.frame(t(vegan_otu(ps_sub)))
  #计算全部均值
  otu_table$mean = rowMeans(otu_table)
  #组合结果
  inde = merge(otu_table,tax_table,by = "row.names", all = TRUE)
  head(inde)
  row.names(inde) = inde$Row.names
  inde$Row.names = NULL

  #添加高级注释信息进入函数
  row.names(inde) = paste("ASV",row.names(inde),sep = "--")
  ##将全部信息合并到节点属性中
  data_plot = merge(vertices_t,inde,by = "row.names",all = TRUE)
  row.names(data_plot) = data_plot$Row.names
  data_plot$Row.names = NULL

  #指定大小映射列将NA值做转换为0
  asa =data_plot$mean
  data_plot$mean[is.na(asa)] = 0
  #选择是否平方根标准化丰度
  # data_plot$mean[!is.na(asa)] = sqrt(data_plot$mean[!is.na(asa)])


  mygraph <- graph_from_data_frame( deg, vertices= data_plot )
  #-----------------------------------设置颜色映射参数-------------------------
  data = create_layout(mygraph, layout = 'circlepack',weight = mean, sort.by = NULL, direction = "out")
  #设置颜色
  data$Phylum = as.character(data$Phylum)
  data$Phylum[is.na(data$Phylum)] = "AA"
  data$Phylum = factor(data$Phylum,levels =unique(data$Phylum) )
  colbar <-length(unique(data$Phylum))
  fil_colors = colorRampPalette(c( "white","#CBD588", "#599861", "orange","#DA5724", "#508578", "#CD9BCD",
                                   "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
                                   "#8569D5", "#5E738F","#D1A33D", "#8A7C64","black"))(colbar)
  names(fil_colors ) = unique(data$Phylum)
  fil_colors[1] = "white"
  p = ggraph(mygraph, layout = 'circlepack',weight = mean, sort.by = NULL, direction = "out") +
    geom_node_circle(aes(fill = as.factor(Phylum),color = as.factor(depth) ) ) +
    scale_fill_manual(values= fil_colors ) +
    scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black", "5"="black", "6"="black", "7"="black") ) +
    geom_node_text( aes(label=new_label), size=6,repel = TRUE) +
    # geom_node_label( aes(label=new_label), size=3,repel = TRUE) +
    theme_void() +
    theme( legend.position="FALSE",plot.margin = unit(rep(0,4), "cm"))#
  # p
  return(list(p,mygraph,fil_colors))
}
