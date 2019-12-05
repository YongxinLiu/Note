rm(list=ls()) 
#导入功能函数
source("./MaptreeForMicrobiome.R")

#------------------otu表格和注释文件都需要矩阵---------------------------------------
#导入otu表格
otu = read.delim("./otutab.txt",row.names = 1)
head(otu)
otu = as.matrix(otu)
# str(otu)
#导入注释文件
tax = read.delim("./taxonomy.txt",row.names = 1)
head(tax)
# tax$ID = row.names(tax)
tax = as.matrix(tax)
#差异表达文件
addtab = read.delim("./KO-WT_all.txt",row.names = 1)

#基本函数，画出maptree
mapdata = data_to_maptree (otu,tax,230)
p1= mapdata[[1]]
p1
ggsave("./maptree1.pdf", p1, width = 12, height =10, device = cairo_pdf, family = "Times New Roman" )

#按照平均丰度修改大小和按照门水平上色颜色
mapadd = maptree_add1_plot(mapdata)
p2 = mapadd[[1]]
#通过gplot修改颜色
p2 +scale_fill_brewer()
# p2 + scale_color_gradient2()

ggsave("./maptree.pdf", p2, width = 12, height =10, device = cairo_pdf, family = "Times New Roman" )



#按照指定列上色
mapadd2 = maptree_add2_plot(mapdata,fillcolor = "level")
#提取作图文件
p3 = mapadd2[[1]]
p3

p3 +scale_fill_brewer()

ggsave("./maptree3.pdf", p3, width = 12, height =10, device = cairo_pdf, family = "Times New Roman" )




