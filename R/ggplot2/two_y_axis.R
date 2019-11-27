# https://blog.csdn.net/enyayang/article/details/98181668

##使用gtable + grid 包中的一些命令，组合ggplot2作图结果，得到双坐标y轴图
library(ggplot2)
library(gtable)
library(grid) 

##定义组合函数
ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  
  # Get the location of the plot panel in g1.
  # These are used later when transformed elements of g2 are put back into g1
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from g2
  index <- which(g2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- g2$grobs[[index]]        # Extract that grob
  ylab <- hinvert_title_grob(ylab)     # Swap margins and fix justifications
  
  # Put the transformed label on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
  index <- which(g2$layout$name == 'axis-l')  # Which grob
  yaxis <- g2$grobs[[index]]          # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(g1)
}

#使用的 R 自带的 mtcars 数据集，分别绘制“mpg-disp”与“mpg-drat”折线图
g1 <- ggplot(mtcars, aes(mpg, disp)) + geom_line(colour = 'blue') + labs(y = 'disp') + theme_bw()
g2 <- ggplot(mtcars, aes(mpg, drat)) + geom_line(colour = 'red') + labs(y = 'drat') + theme_bw() %+replace% 
  theme(panel.background = element_rect(fill = NA))

#使用 ggplot2.two_y_axis() 整合结果，并保存图片
ggsave('g1.pdf', g1)
ggsave('g2.pdf', g2)
ggsave('g12.pdf', ggplot2.two_y_axis(g1, g2))
