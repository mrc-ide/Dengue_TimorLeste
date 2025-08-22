#function to plot observed (red) and estimated (green seroprevalence)
plot_f <- function(z_statA, true_IC_serop, EA_name) {
  
  plot <- ggplot() +
    
    geom_line(data=z_statA, aes(x=z_statA$age, y=z_statA$med*100), 
              color="black")+
    
    geom_ribbon(data=z_statA, aes(x= z_statA$age, y=z_statA$med*100, 
                                  ymin = z_statA$ciL*100, ymax = z_statA$ciU*100), 
                color ="green", fill="green", alpha =0.5)+
    
    geom_point(data = true_IC_serop, aes(x=true_IC_serop$age, y=true_IC_serop$med*100),
               size = 2, color="red")+
    
    geom_errorbar(data = true_IC_serop, aes(x=true_IC_serop$age, y=true_IC_serop$med*100,
                                            ymin = true_IC_serop$ciL*100, ymax = true_IC_serop$ciU*100),
                  color="red", width= 0.5)+
    
    ylim(0,100) +
    theme_bw()+
    labs(x = element_blank(), y = element_blank())+
    theme(plot.title = element_text(hjust = 0.5))+ 
    # ggtitle(EA_name)+
    facet_wrap(~ EA, scales = "free") +
    theme(legend.position = "none")+
    xlab("age") + ylab("% seroprevalence")
  
  return(plot)
}