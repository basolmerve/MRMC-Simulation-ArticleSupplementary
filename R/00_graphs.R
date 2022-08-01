library(dplyr)
library(ggplot2)
library(magrittr)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


veri1 <- read.table("veri.txt", header = TRUE)
str(veri1)

# veri1$MU <- as.factor(veri1$MU)
# veri1$READER <- as.factor(veri1$READER)
# veri1$N <- as.factor(veri1$N)

# Filtering
nR <- 10
varEqual <- 2  # 1: Equal, 2: Unequal

if (varEqual == 2){
  rmar <- 0
} else {
  rmar <- 65
}

veri <- filter(veri1, VAR_EQUAL == varEqual, READER == nR, Model != "BootT")

plotTitle <- paste(ifelse(varEqual == 1, " Equal", " Unequal"), " Variance Components", sep = "")
plotData <- transform(veri, 
                      N = factor(N, levels = c(25, 50, 100),
                                 labels = c("n = 25", "n = 50", "n = 100")),
                      MU = factor(veri$MU, levels = c(2.5, 1.5, 0.75), 
                                  labels = c("AUC = 0.961", "AUC = 0.856", "AUC = 0.702")),
                      VAR_EQUAL = factor(veri$VAR_EQUAL, levels = c(1, 2), 
                                         labels = c("Equal", "Unequal")))

plotData$Model <- factor(plotData$Model, 
                         levels = c("DBM", "OR", "MM", "BCa", "Percentile"))

p2 <- ggplot(plotData, aes(y = TypeI_Err, x = CorStr, shape = Model, color = Model)) + 
  geom_hline(yintercept = c(0.06), size = 0.2, colour = "gray40", linetype = 2) +
  geom_hline(yintercept = c(0.05), size = 0.2, linetype = 3) +
  geom_hline(yintercept = c(0.04), size = 0.2, colour = "gray40", linetype = 2) +
  geom_point(size = 1.5) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(face = "plain", margin = margin(0, 0, 5, 0), hjust = 0.5),
        axis.title.x = element_text(margin = margin(5,0,5,0), face = "plain"),
        axis.title.y = element_text(margin = margin(0,5,0,5), face = "plain"),
        axis.text.x = element_text(margin = margin(2,0,0,0), face = "plain"),
        axis.text.y = element_text(margin = margin(0,2,0,0), face = "plain"), 
        plot.margin = margin(0,rmar,0,0)) +
  xlab("Correlation structure") +
  ylab("Type-I error rate") +
  ggtitle(plotTitle) + 
  ylim(c(0, 0.1)) +
  facet_grid(MU ~ N) + 
  scale_shape_manual(values = c(0,1,2,3,5)) + 
  scale_color_manual(values = c("gray25","gray25","gray25","gray50","gray50"))

if (varEqual == 1){
  p1 <- p1 + guides(shape = FALSE, color = FALSE)
}

cairo_pdf(file = "figures/10Readers.pdf", height = 4, width = 10, pointsize = 6)
multiplot(p1, p2, cols = 2)
dev.off()

ggsave(plot = multiplot(p1, p2, cols = 2), 
       filename = "figures/10Readers.eps", device = "eps", width = 10, height = 4, units = "in", pointsize = 6)

png(filename = "figures/10Readers.png", width = 3000, height = 1200, res = 300, pointsize = 6)
multiplot(p1, p2, cols = 2)
dev.off()

#### Confidence Intervals
ciData <- read.table("confInts.txt", header = TRUE)
str(ciData)

ciData <- mutate(ciData, mid = (Lower + Upper) / 2)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.5) # move them .05 to the left and right

ggplot(ciData, aes(x = Method, y = mid, colour = Model)) + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, position = pd, size = 0.8) + 
  ylim(c(-0.1, 0.05)) + 
  theme_bw(base_size = 14) + 
  theme(panel.grid = element_blank(), plot.title = element_text(face = "plain", margin = margin(5, 0, 10, 0), hjust = 0.5),
        axis.title.x = element_text(margin = margin(15,0,10,0), face = "plain"),
        axis.title.y = element_text(margin = margin(0,15,0,10), face = "plain")) + 
  geom_hline(yintercept = c(0), size = 0.2, linetype = 3) + 
  geom_point(aes(y = Lower), position = pd, size = 2) + 
  geom_point(aes(y = Upper), position = pd, size = 2) +
  xlab("AUC Estimation Method") + 
  ylab(expression("95% CI for "("AUC"[1] - "AUC"[2])))
# scale_colour_manual(values = rep("black", 5))



## Only Fixed:
plotData <- filter(ciData, Method == 'Fixed')

pd <- position_dodge(0) # move them .05 to the left and right

ggplot(plotData, aes(x = Model, y = mid, colour = Model)) + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, size = 0.8) + 
  ylim(c(-0.1, 0.05)) + 
  theme_bw(base_size = 14) + 
  theme(panel.grid = element_blank(), plot.title = element_text(face = "plain", margin = margin(5, 0, 10, 0), hjust = 0.5),
        axis.title.x = element_text(margin = margin(15,0,10,0), face = "plain"),
        axis.title.y = element_text(margin = margin(0,15,0,10), face = "plain")) + 
  geom_hline(yintercept = c(0), size = 0.2, linetype = 3) + 
  geom_point(aes(y = Lower), size = 2) + 
  geom_point(aes(y = Upper), size = 2) +
  xlab("AUC Estimation Method") + 
  ylab(expression("95% CI for "("AUC"[1] - "AUC"[2])))

