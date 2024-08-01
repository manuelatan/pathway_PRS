library(dplyr)
library(ggplot2)

#Function to make forest plots for pathway-PRS with all the pathways
#x is the dataframe
#Does not include meta-analysis results
makeForestPlot <- function(x) {
  
  #Plot forest plot
  
  #use the levels = rev() to reverse order of pathway factors
  #Otherwise they are shown in alphabetical order from bottom to top
  #We want alphabetical order top to bottom
  plot <- ggplot(data = x, aes(x = RE_SMD, y = factor(pathway, 
                                                      levels = rev(levels(factor(pathway)))))) +
    
    #Shape, size and colour of points (betas)
    geom_point(shape = 16, colour = "blue", size = 3) +
    
    #Error bars - 95% confidence interval
    geom_errorbarh(aes(xmin = X95CI_lower, xmax = X95CI_upper), height = 0.0, colour = "blue") +
    
    #Dashed line through beta = 0
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    #Axis labels
    xlab("Beta coefficient (95% Confidence Interval)") +
    ylab("Pathway") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = "bold"),
          plot.title = element_text(size = 18))
  
  print(plot)
  
  
}



#Function to make forest plots for pathway-PRS for a specific pathway
#x is the dataframe
#Plots cohort results separately for each pathway and outcome, and plots meta-analysis results
makeForestPlot_cohorts <- function(x) {
  
  #Plot forest plot
  x <- x %>% 
    mutate(meta = ifelse(dataset == "Random effects meta-analysis", 1, 0))
  #use the levels = rev() to reverse order of pathway factors
  #Otherwise they are shown in alphabetical order from bottom to top
  #We want alphabetical order top to bottom
  #Random effects meta-analysis is at the bottom
  reshuffled <- fct_relevel(x$dataset, "Random effects meta-analysis", after = Inf)
  
  
  plot <- ggplot(data = x, aes(x = Coeff, y = factor(dataset, 
                                                        levels = rev(levels(reshuffled))))) +
    
    #Shape, size and colour of points (betas)
    geom_point(aes(shape = factor(meta), fill = factor(meta), colour = factor(meta)), size = 3) +
    scale_shape_manual(values=c(21,23)) +
    scale_color_manual(values=c("blue", "red")) +
    scale_fill_manual(values = c("blue", "red")) +
    
    #Error bars - 95% confidence interval
    geom_errorbarh(aes(xmin = (Coeff - 1.96*se), xmax = (Coeff + 1.96*se), colour = factor(meta)), height = 0.0, ) +
    
    #Dashed line through beta = 0
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    #Axis labels
    xlab("Beta coefficient (95% Confidence Interval)") +
    ylab("Cohort") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = "bold"),
          legend.position = "none")
  
  print(plot)
  
  
}




#Function to make forest plots for progression PRS, separated by cohort
#x is the dataframe
#Plots cohort results separately for each pathway and outcome, and plots meta-analysis results
makeProgressionPRSplot <- function(x) {
  
  #Plot forest plot
  x <- x %>% 
    mutate(meta = ifelse(dataset == "Random effects meta-analysis", 1, 0))
  #use the levels = rev() to reverse order of pathway factors
  #Otherwise they are shown in alphabetical order from bottom to top
  #We want alphabetical order top to bottom
  #Random effects meta-analysis is at the bottom
  reshuffled <- fct_relevel(x$dataset, "Random effects meta-analysis", after = Inf)
  
  
  plot <- ggplot(data = x, aes(x = Coeff, y = factor(dataset, 
                                                     levels = rev(levels(reshuffled))))) +
    
    #Shape, size and colour of points (betas)
    geom_point(aes(shape = factor(meta), fill = factor(meta), colour = factor(meta)), size = 3) +
    scale_shape_manual(values=c(21,23)) +
    scale_color_manual(values=c("blue", "red")) +
    scale_fill_manual(values = c("blue", "red")) +
    
    #Error bars - 95% confidence interval
    geom_errorbarh(aes(xmin = (Coeff - 1.96*se), xmax = (Coeff + 1.96*se), colour = factor(meta)), height = 0.0, ) +
    
    #Dashed line through beta = 0
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    #Axis labels
    xlab("Beta coefficient (95% Confidence Interval)") +
    ylab("Cohort") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = "bold"),
          legend.position = "none")
  
  print(plot)
  
  
}


#Function to make forest plots for progression PRS, without meta-analysis results separated by cohort
#x is the dataframe
makeProgressionPRSplot_nometa <- function(x) {
  
  plot <- ggplot(data = x, aes(x = Coeff, y = factor(dataset))) +
    
    #Shape, size and colour of points (betas)
    geom_point(shape = 21, fill = "blue", colour = "blue", size = 3) +
    
    #Error bars - 95% confidence interval
    geom_errorbarh(aes(xmin = (Coeff - 1.96*se), xmax = (Coeff + 1.96*se)), colour = "blue", height = 0.0) +
    
    #Dashed line through beta = 0
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    #Axis labels
    xlab("Beta coefficient (95% Confidence Interval)") +
    ylab("Cohort") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = "bold"),
          legend.position = "none") +
    
    #Set x axis limits
    xlim(-0.1, 0.6)
  
  print(plot)
  
  
}

