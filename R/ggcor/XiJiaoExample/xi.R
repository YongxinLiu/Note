library(ggplot2)
spec <- readr::read_csv(file = file.choose())
env <- readr::read_csv(file = file.choose())
corr <- fortify_cor(env, type = "upper", cor.test = TRUE)
corr2 <- fortify_cor(spec, env, cor.test = TRUE) %>% 
  filter(p.value < 0.05) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0, Inf),
                  labels = c("< 0", ">= 0"),
                  right = FALSE))

quickcor(corr) +
  geom_square(colour = NA) +
  geom_cross(colour = "grey60", size = 0.4) +
  anno_link(aes(x = x + 1, colour = rd, linetype = rd, size = abs(r)), data = corr2,
            node.size = c(6, 5), node.shape = c(24, 21),
            node.colour = NA, node.fill = c("red", "grey60")) +
  anno_link_label(size = 4.5, nudge_x = -0.5) +
  scale_fill_gradient2n(colours = rev(red_blue()), 
                        breaks = c(-1, -0.5, 0, 0.5, 1),
                        limits = c(-1, 1)) +
  scale_color_manual(values = c("#FF7F00", "#4DAF4A")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_size(range = c(0.25, 2), breaks = c(0.4, 0.5, 0.6, 0.7, 0.8), limits = c(0.4, 0.8)) +
  guides(size = guide_legend(override.aes = list(colour = "grey35"), order = 3),
         fill = guide_colorbar(title = "Corr", order = 100))
ggsave("xi.pdf", width = 13, height = 9, dpi = 600)
