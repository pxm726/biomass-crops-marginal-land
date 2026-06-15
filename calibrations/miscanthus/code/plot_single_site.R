library(ggplot2)
Compare_Observed_and_Predicted <- function(observeddata, predicteddata, xlabtitle, ylabtitle){
  
  # Create a single factor variable for legend grouping
  observeddata$varname <- as.factor(observeddata$varname)
  predicteddata$varname <- as.factor(predicteddata$varname)
  
  compareplot <- ggplot(data = observeddata)
  compareplot <- compareplot + geom_point(aes(x = x, y = y, color = varname), size = 3)
  compareplot <- compareplot + geom_line(data = predicteddata, aes(x = x, y = y, color = varname), linewidth = 1)
  
  compareplot <- compareplot + xlab(xlabtitle)
  compareplot <- compareplot + ylab(ylabtitle)
  
  N <- length(unique(observeddata$varname))
  
  if(N == 1){
    compareplot <- compareplot + scale_colour_manual(values = c("black"))
  } else {
    # Professional color palette that works well on white background
    cbPalette <- c("#117733", "#999933", "#882255", "#E69F00", "#009E73", "#332288", "#0072B2", "#D55E00")
    cbPalette <- cbPalette[1:N]
    compareplot <- compareplot + scale_colour_manual(values = cbPalette)
  }
  
  # Single legend with both points and lines
  compareplot <- compareplot + 
    guides(colour = guide_legend(
      title = NULL,
      override.aes = list(
        shape = 16,        # Points for all
        linetype = 1,      # Lines for all
        linewidth = 1, 
        size = 3
      )
    ))
  
  # Apply publication-ready theme
  compareplot <- compareplot +
    theme_classic(base_size = 12) +
    theme(
      # White background
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Grid lines (subtle, optional - remove if not desired)
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      
      # Axes
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 13, face = "bold", color = "black"),
      
      # Legend
      legend.background = element_rect(fill = "white", color = "NA", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "bottom",  # Inside plot area
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(0.6, "cm"),
      
      # Plot margins
      plot.margin = margin(15, 15, 15, 15)
    )
  
  return(compareplot)  
}


# library(ggplot2)
# Compare_Observed_and_Predicted <- function(observeddata,predicteddata,xlabtitle,ylabtitle){
#   
#   observeddata$pointlegend <- as.factor(paste("observed",observeddata$varname,sep =" "))
#   predicteddata$linelegend <- as.factor(paste("predicted", predicteddata$varname, sep = " "))
#   
#   compareplot <- ggplot(data=observeddata)
#   compareplot <- compareplot + geom_point(aes(x=x,y=y,color=pointlegend))
#   compareplot <- compareplot + geom_line(data = predicteddata, aes(x=x,y=y,color=linelegend))
#   
#   
#   compareplot <- compareplot + xlab(xlabtitle)
#   compareplot <- compareplot + ylab(ylabtitle)
#   
#   N <- length(unique(observeddata$varname))
#   if(N==1){
#     compareplot <- compareplot + scale_colour_manual(values= c("black","black"))
#     shapelist <- c(16,NA)
#     linelist <- c(0,1)
#   } else {
#     cbPalette <- c("#000000", "#E69F00", "#CC79A7", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
#     cbPalette <- cbPalette[1:N]
#     compareplot <- compareplot + scale_colour_manual(values = rep(cbPalette,2))
#     shapelist <- c(rep(16,N),rep(NA,N))
#     linelist <- c(rep(0,N),rep(1,N))
#   }
#   compareplot <- compareplot + guides(colour = guide_legend(title=NULL,override.aes = list(shape=shapelist,linetype=linelist)),shape=guide_legend(title=NULL)) 
#   return(compareplot)  
# }
