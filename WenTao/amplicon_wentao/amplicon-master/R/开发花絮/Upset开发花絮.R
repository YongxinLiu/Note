

# 检查R包安装关系与依赖
# 依赖关系检测与安装
p_list = c("ggplot2", "dplyr", "multcompView") # "agricolae"
for(p in p_list){
  if (!requireNamespace(p)){
    install.packages(p)}
  suppressPackageStartupMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
}

# 
# library("rjava")
# 
# library(BiocManager)
# install("rjava")
require(ggplotify)


### 默认版本的upset
p <- upset(ven2, sets = colnames(ven2),
          number.angles = 30, point.size = 2, line.size = 1,
          mainbar.y.label = "OTU", sets.x.label = "OTU Per Treatment",
          text.scale = c(2, 2, 2,2, 2, 2),mb.ratio = c(0.7, 0.3),order.by = "freq",keep.order = TRUE,
          queries = list(list(query = intersects, params = 
                                list(colnames(ven2)), color = "red", active = T),
                         list(query = intersects, params = 
                                list(colnames(ven2)), color = "red", active = T),
                         list(query = intersects, params = 
                                list(colnames(ven2)), color = "red", active = T)))
# p1 <- as.ggplot(p)



p
library("gridExtra")
gridExtra:::grid.newpage(p)


library(tidyverse)
library(ggupset)

tidy_movies %>%
  distinct(title, year, length, .keep_all=TRUE) %>%
  ggplot(aes(x=Genres)) +
  geom_bar() +
  scale_x_upset(n_intersections = 20)


aa =as.data.frame(tidy_movies)
head(aa)

aa$Genres
bb = distinct(aa,title, year, length, .keep_all=TRUE)
head(bb)
bb$Genres


library(tibble)
library(tidyr)
library(dplyr)

aa = t(ven2)
head(aa)
## 调整数据格式
tidydata <- aa %>%
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  gather(OTU, Member, -ID) %>%
  filter(Member==1) %>%
  select(- Member) %>% 
  group_by(OTU) %>%
  summarize(ID = list(ID))

ggplot(tidydata,aes(x = ID)) +
  geom_bar() +
  scale_x_upset()


#--------------测试 案例
library(UpSetR)

movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=TRUE, sep=";" )

require(ggplot2); require(plyr); require(gridExtra); require(grid);

between <- function(row, min, max){
  newData <- (row["ReleaseDate"] < max) & (row["ReleaseDate"] > min)
}

plot1 <- function(mydata, x){
  myplot <- (ggplot(mydata, aes_string(x= x, fill = "color"))
             + geom_histogram() + scale_fill_identity()
             + theme(plot.margin = unit(c(0,0,0,0), "cm")))
}

plot2 <- function(mydata, x, y){
  myplot <- (ggplot(data = mydata, aes_string(x=x, y=y, colour = "color"), alpha = 0.5)
             + geom_point() + scale_color_identity()
             + theme_bw() + theme(plot.margin = unit(c(0,0,0,0), "cm")))
}

attributeplots <- list(gridrows = 55,
                       plots = list(list(plot = plot1, x= "ReleaseDate",  queries = FALSE),
                                    list(plot = plot1, x= "ReleaseDate", queries = TRUE),
                                    list(plot = plot2, x = "ReleaseDate", y = "AvgRating", queries = FALSE),
                                    list(plot = plot2, x = "ReleaseDate", y = "AvgRating", queries = TRUE)),
                       ncols = 3)

upset(movies, nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))

upset(movies, sets = c("Drama", "Comedy", "Action", "Thriller", "Western", "Documentary"),
      queries = list(list(query = intersects, params = list("Drama", "Action")),
                     list(query = between, params = list(1970, 1980), color = "red", active = TRUE)))

upset(movies, attribute.plots = attributeplots,
      queries = list(list(query = between, params = list(1920, 1940)),
                     list(query = intersects, params = list("Drama"), color= "red"),
                     list(query = elements, params = list("ReleaseDate", 1990, 1991, 1992))),
      main.bar.color = "yellow")





