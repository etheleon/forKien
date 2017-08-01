R Code required to replicate the following plot

![thesis](https://github.com/etheleon/forKien/blob/master/wesleyPlot.png)

```{r}
library(tidyverse)
library(magrittr)
library(cowplot)
load("data.rda")

cutoff = 0.4

koDF2 %<>% filter(cum < cutoff)
locDF %<>% filter(cum < cutoff)
rect = locDF[1,]
rect %<>% mutate(y= 0.5, ymin = unique(locDF$msaS), ymax = unique(locDF$msaE), xmin = 0, xmax=cutoff, type="loc")
koi = 'K00927'
nameOfko = 'phosphoglycerate kinase [EC:2.7.2.3]'

combinedPlot = ggplot(data=dummyDataFrame, mapping = aes(x=cum, y=y)) +
    facet_grid(type~., labeller = labeller(type=c(expression="Expression (Relative percentage)", loc="MSA Coordinates")), scale="free", switch="y") +
    geom_line(data=koDF2, aes(x=cum, y=count.gDNA.perc), size=1) +
    geom_segment(data=koDF2, aes(x=cum, xend = cum, y=count.gDNA.perc, yend=count.gDNA.perc+count.cDNA.perc, color=status, alpha=spanning))+
    geom_point(data=koDF2, aes(x=cum, y=count.gDNA.perc, size=sqrt(count.cDNA.perc), color =status, alpha=spanning)) +
    geom_rect(data=rect, aes(ymin = ymin, ymax = ymax, xmin = 0, xmax=cutoff), alpha=0.3) +
    geom_segment(data=locDF, aes(x=cum, xend=cum, y=start, yend=end, color=status, alpha=spanning)) +
    scale_x_continuous(labels=scales::percent)+
    scale_color_manual("Taxonomy (LCA)",
        breaks = c("genus", "higher", "unclassifiable"),
        labels=c("Genus Rank", "Above Genus", "Unclassifiable"),
        drop=FALSE, values=c("#aac476", "#d16461", "#5b95c4"))+
    scale_size_continuous("Expression", labels= paste0(as.character((seq(0, 0.20, 0.05)^2)*100), "%"),breaks=seq(0, 0.20, 0.05))+
    scale_alpha_continuous("Spanning MDR", breaks=c(0.35, 1), labels=c("Not spanning", "Spanning")) +
    ggtitle(sprintf("%s - %s", koi, nameOfko))+
    theme(strip.background = element_blank()) +
    ylab(NULL) +
    xlab("Cumulative abundance of top OTUs in MDR region")

yticks = ggplot_build(combinedPlot)$layout$panel_ranges[[1]]$y.major_source
yticksSmall = ggplot_build(combinedPlot)$layout$panel_ranges[[2]]$y.major_source
yticks = sort(c(yticks, as.integer(unique(locDF$msaS)), as.integer(unique(locDF$msaE)), yticksSmall))

ylabels = sapply(unique(yticks), function(x) {
    if(x==0) { "0"}
    else if(x < 1) { paste0(as.character(x * 100), "%")}
    else{as.character(as.integer(x))}
})
combinedPlot = combinedPlot + scale_y_continuous(breaks=yticks, labels=ylabels)

shannonPlot = shannon %>% filter(ko == koi) %>% filter(loc >= start, loc <= end) %>%
    ggplot(aes(x=loc, y=`shannon Score`, group=ko)) + geom_line() +
    theme(panel.background = element_rect(fill = "#D3D3D350")) + ylab("Shannon Entropy") + xlab("MSA Coordinates")

xbreaks = ggplot_build(shannonPlot)$layout$panel_ranges[[1]]$x.major_source
xbreaks = xbreaks[-1]
xbreaks = xbreaks[-length(xbreaks)]
xbreaks = c(xbreaks, as.integer(unique(locDF$msaS)), as.integer(unique(locDF$msaE)))
shannonPlot = shannonPlot + scale_x_continuous(breaks = xbreaks)

combinedPlot = ggdraw() +
    draw_plot(combinedPlot, 0, 0, 1, 1) +
    geom_rect(xmin = 0.54, xmax = 0.85, ymin = 0.07, ymax =  0.28, fill = "white", color="black") +
    draw_plot(shannonPlot, 0.55, 0.07, 0.3, 0.2) +
    draw_plot_label(c("A", "B"), c(0, 0.54), c(1, 0.28), size = 15)

pdf("wesleyPlot.pdf", w=15, h=10)
combinedPlot
dev.off()
ggsave(combinedPlot, file=”wesleyPlot.png”, w=15, h=10)
```
