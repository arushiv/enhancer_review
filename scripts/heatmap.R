library(ggplot2)
library(reshape2)
library(scales)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=T)

#EnhancerType1		 Cell_type1	EnhancerType2	Cell_type2	toal_length1	base_overlap	total_regions1	region_overlap

d$base_fraction <- d$base_overlap/d$total_length1
d$region_fraction <- d$region_overlap/d$total_regions1

levels(d$EnhancerType1)
levels(d$EnhancerType2)

d$EnhancerType1 <- factor(d$EnhancerType1, levels=c("Typical_Enhancers","Super_Enhancers","Stretch_Enhancers"))
levels(d$EnhancerType1)

makePlot <- function(d){
    p <- ggplot(d, aes(x = Cell_type2, y =Cell_type1)) +
        geom_tile(aes(fill=region_fraction)) +
        facet_grid(EnhancerType1~EnhancerType2) +
        scale_fill_gradientn(colours=c("white","yellow", "orange", "red"),
                             values=rescale(c(0,.2,.5,.75)),
                             name="Fraction of regional overlap") +
        theme(text = element_text(size=10),
              axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=10),
              axis.text.y = element_text(size=10),
              panel.background = element_rect(fill="white"))

    return(p)
}

pdf(args[2], height=7, width=10)
makePlot(d)
dev.off()


