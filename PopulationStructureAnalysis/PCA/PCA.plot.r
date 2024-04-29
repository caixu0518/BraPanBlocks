library(ggplot2)

data <- read.table("PCA.results.format.lst", header=T, sep="\t")

custom_colors <- c("Oilseed" = "#6e6b41", "Grelos" = "#000000", "Turnip" = "#278EAA", "JPNLeafy" = "#FFC125", "Taicai" = "#DC7017", "Caixin" = "#40de5a", "Caixin1" = "#40de5a", "Caixin2" = "#40de5a", "Pakchoi" = "#16a951", "Zicaitai" = "#B08ABD", "Wucai" = "#9ed900", "Wutacai" = "#057748", "ChineseCabbage" = "#ff2d51", "Rapini" = "#700961")

p1 <- ggplot(data, aes(x=pca1, y=pca2, fill=Group, color=Group)) + 
      geom_point(size=1.6, alpha=0.5, stroke=0) +
      scale_fill_manual(values = custom_colors) +
      scale_color_manual(values = custom_colors) +
      theme(panel.background = element_blank(),axis.line = element_line(), axis.text.x=element_text(size=6, angle=0), legend.key.size = unit(0.4,'cm'), legend.text = element_text(size=6,angle=0)) +
      xlab("PC1") + ylab("PC2")

ggsave(file="PCA.dotplot.pdf", plot=p1, width = 5.5, height = 4)
