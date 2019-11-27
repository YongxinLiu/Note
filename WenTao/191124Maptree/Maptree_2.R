
#------------------otu表格和注释文件都需要矩阵---------------------------------------
#导入otu表格
otu = read.delim("./otutab.txt",row.names = 1)
head(otu)
otu = as.matrix(otu)
str(otu)
#导入注释文件
tax = read.delim("./taxonomy.txt",row.names = 1)
head(tax)
tax = as.matrix(tax)
# taxa_names(tax)
#导入分组文件
map = read.delim("./metadata.tsv",row.names = 1)
head(map)

#差异表格添加
KOWT = read.delim("./KO-WT_all.txt",row.names = 1)

# 防止合并后名称改变
colnames(KOWT) = paste("KOWT",colnames(KOWT),sep = "_")
head(KOWT)

N = 100
#选择需要的列用于上色
fillcolor = "KOWT_level"

##-------------------------------------------依赖包-----------
library(ggraph)
library(igraph)
library(tidyverse)
library(viridis)
library(data.tree)
library(phyloseq)
library(ggtree)
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
               sample_data(map) ,
               tax_table(tax)
               
)
ps

# 命名并保存
# name = "liu"
# path = "./"
# proID = paste(path,"ps_",name,".rds",sep = "")

# saveRDS(ps,proID)
ps1_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela 
# 提取otu表格

otu_table = as.data.frame(t(vegan_otu(ps1_rela)))
otu_table$mean = rowMeans(otu_table)
otu_table$ID = row.names(otu_table)
head(otu_table)
#按照从大到小排序

otu_table<- arrange(otu_table, desc(mean))
subtab = head(otu_table,N)
dim(subtab)
# row.names(subtab) = subtab$ID
#对phyloseq取子集
ps_sub <- ps1_rela %>%
  subset_taxa(
    
    row.names(tax_table(ps1_rela ))%in%subtab$ID
  )
ps_sub


#--------------------构造边和节点文件-----------------
#注意我们本次制作的是OTU水平的maptree

tax = as.data.frame(vegan_tax(ps_sub))
tax_point = as.data.frame(vegan_tax(ps_sub))
head(tax)
#添加水平识别标示
tax$Kingdom = paste("D_",tax$Kingdom,sep = "--")
tax$Phylum = paste(tax$Kingdom,"P_",tax$Phylum,sep = "--")
tax$Class = paste(tax$Phylum,"C_",tax$Class,sep = "--")
tax$Order = paste(tax$Class,"O_",tax$Order,sep = "--")
tax$Family = paste(tax$Order,"F_",tax$Family ,sep = "--")
tax$Genus = paste(tax$Family,"G_",tax$Genus,sep = "--")
tax$Species = paste(tax$Genus,"S_",tax$Species,sep = "--")
# str(tax)#保证是字符型变量
# # 对不同分类水平的tax未知物种进行区分
#设置未识别物种标示

#开始构造分组文件
tax_table(ps_sub) = as.matrix(tax)
tax$OTU = paste("ASV",row.names(tax),sep = "--")
#提取分组信息
Desep_group <- as.character(colnames(tax))
Desep_group
head(tax)
#构造两个分组的全部组合
# aaa = combn(Desep_group,2)
#按照form--to格式来排布数据
i = 1
edge_t = tax[c(Desep_group[1:2])]
dim(edge_t)
edge_t <-  dplyr::distinct(edge_t)
dim(edge_t)
i= 7
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
  dim(deg)
  
}
deg = edge_t


tree <- FromDataFrameNetwork(deg)

tree
vertices_t  <-  data.frame(
  name = unique(c(as.character(deg$from), as.character(deg$to)))
)

head(vertices_t )
# Then I can easily get the level of each node, and add it to the initial data frame:
mylevels <- data.frame( name=tree$Get('name'), level=tree$Get("level") )
vertices_t <- vertices_t %>% 
  left_join(., mylevels, by=c("name"="name"))
dim(vertices_t)


# sapply(strsplit(basename(vertices_t$name), "_"), `[`, 3)

AA = rep("A",length(vertices_t$name))
for (i in 1:length(vertices_t$name)) {
  AA [i] = strsplit(basename(as.character(vertices_t$name)), "--")[[i]][length(strsplit(basename(as.character(vertices_t$name)), "--")[[i]])]
  
}

vertices_t$shortName<- AA
# #提取tax前标示作为level
# vertices_t$level = sapply(strsplit(basename(vertices_t$name), "_"), `[`, 1)
# 

vertices_t$level



vertices_t <- vertices_t %>%
  mutate(new_label=ifelse(level==2, shortName, NA))



##下面我想将不同节点的丰度信息纳入size指标。

tax_table = as.data.frame(vegan_tax(ps_sub))

sub_design =as.data.frame( sample_data(ps_sub))
otu_table = as.data.frame(t(vegan_otu(ps_sub)))
# head(otu_table)
count = otu_table
count2 = as.data.frame(t(count))
head(count2)
library("tibble")
#数据分组
iris.split <- split(count2,as.factor(sub_design$Group))
#数据分组计算平均值
iris.apply <- lapply(iris.split,function(x)colMeans(x[]))
# 组合结果
iris.combine <- do.call(rbind,iris.apply)
otu_mean = t(iris.combine)
# head(otu_mean)
colnames(otu_mean) = paste("mean_",colnames(otu_mean),sep = "")

#计算全部均值
otu_table$mean = rowMeans(otu_table)
#组合结果
index = merge(otu_table,tax_table,by = "row.names", all = TRUE)
# head(index)
row.names(index) = index$Row.names
index$Row.names = NULL
#合并分组均值文件
inde =  merge(index,otu_mean,by = "row.names", all = TRUE)
head(inde)

row.names(inde) = inde$Row.names
inde$Row.names = NULL
head(inde)

dim(inde)
inde = merge(inde,KOWT,by = "row.names",all.x = TRUE)
row.names(inde) = inde$Row.names
inde$Row.names = NULL
head(inde)



colnames(inde) <- gsub(fillcolor,"fill_color",colnames(inde))


row.names(inde) = paste("ASV",row.names(inde),sep = "--")
# write.csv(inde,"./inde.csv")

head(vertices_t)
row.names(vertices_t) = vertices_t$name

##将全部信息合并到节点属性中
data_plot = merge(vertices_t,inde,by = "row.names",all = TRUE)
row.names(data_plot) = data_plot$Row.names
data_plot$Row.names = NULL
dim(data_plot )


#开始设置大小
asa =data_plot$mean
data_plot$mean[is.na(asa)] = 0
#是否标准化丰度
# data_plot$mean[!is.na(asa)] = sqrt(data_plot$mean[!is.na(asa)])

head(data_plot)

data_plot$KOWT_level

mygraph <- graph_from_data_frame( deg, vertices= data_plot )


data = create_layout(mygraph, layout = 'circlepack',weight = "mean", sort.by = NULL, direction = "out")

#设置颜色

i = "fill_color"


head(data)


data[[i]] = as.character(data[[i]])
data[[i]][is.na(data[[i]])] = "AA"
unique(data[[i]])
data[[i]] = factor(data[[i]],levels =unique(data[[i]]) )
levels(data[[i]])

data[i] = data[[i]]



library(RColorBrewer)#调色板调用包
mi = display.brewer.pal(8,"Set1")
colbar <-length(unique(data[[i]]))
library("scales")
# fil_colors = colorRampPalette(c( "white","#CBD588", "#599861", "orange","#DA5724", "#508578", "#CD9BCD",
#                                  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
#                                  "#8569D5", "#5E738F","#D1A33D", "#8A7C64","black"))(colbar)
mi = brewer.pal(9,"Set1")
fil_colors = mi[1:colbar ]
names(fil_colors ) = unique(data$KOWT_level)
library("scales")
# show_col(fil_colors)
fil_colors

write.csv(data,"data2.csv")
p = ggraph(mygraph, layout = 'circlepack',weight = "mean", sort.by = NULL, direction = "out") + 
  geom_node_circle(aes(fill = as.factor(fill_color),color = as.factor(depth) ) ) +
  scale_fill_manual(values= fil_colors ) +
  scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black", "5"="black", "6"="black", "7"="black") ) +
  geom_node_text( aes(label=new_label), size=6,repel = TRUE) +
  # geom_node_label( aes(label=new_label), size=3,repel = TRUE) +
  theme_void() + 
  theme( legend.position="FALSE",plot.margin = unit(rep(0,4), "cm"))#
p

ggsave("./maptree3.pdf", p, width = 12, height =10, device = cairo_pdf, family = "Times New Roman" )


p = ggraph(mygraph, layout = 'circlepack',weight = "mean", sort.by = NULL, direction = "out") + 
  geom_node_circle(aes(fill = as.factor(fill_color),color = as.factor(depth) ) ) +
  scale_fill_manual(values= fil_colors ) +
  scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black", "5"="black", "6"="black", "7"="black") ) +
  # geom_node_text( aes(label=new_label), size=3,repel = TRUE) +
  geom_node_label( aes(label=new_label), size=6,repel = TRUE) +
  theme_void() + 
  theme( legend.position="FALSE",plot.margin = unit(rep(0,4), "cm"))#
p

ggsave("./maptree4.pdf", p, width = 12, height =10, device = cairo_pdf, family = "Times New Roman" )






