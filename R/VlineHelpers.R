# previous one ----
PlotVline <- function(all_nonB, metric, nonb_type, x_title){

  mr <- all_nonB %>% filter(motif_type == nonb_type)

  z = mr %>% group_split(group) %>%
    lapply(function(tib){pull(tib,metric)})
  z0=z[[1]]
  z1=z[[2]]

  ggplot()+
    geom_vline(xintercept = z0, color="grey80")+
    geom_density(aes(z0), color="black")+
    #geom_vline(xintercept = qnorm(c(0.05, 0.95),mean(z0),sd(z0)), color="blue")+
    geom_vline(xintercept = z1,color="red", size =0.5)+
    #ggprism::theme_prism()+
    theme_linedraw()+
    theme(axis.text  = element_text(size = 10, face="bold"),
          axis.title  = element_text(size = 10, face="bold"),
          plot.title = element_text(size=12, face = "bold")) +
    theme(panel.grid.major = element_line(linetype = "blank"),
          panel.grid.minor = element_line(linetype = "blank"))+
    labs(title= nonb_type, x=x_title)
}
PlotVline_clean <- function(all_nonB,
                            metric="total_overlapped_region",
                            nonb_type = "Z_DNA_Motif",
                            x_title=NULL){

  mr <- all_nonB %>% filter(motif_type == nonb_type)

  z = mr %>% group_split(group) %>%
    lapply(function(tib){pull(tib,metric)})
  z0=z[[1]]
  z1=z[[2]]

  p_value <- mean(z0 > z1)

  if (p_value > 0.05) {
    vline_col <- "blue"
    vline_size = 0.5
    vline_type = "dashed"
  } else  {
    vline_col <- "red"
    vline_size <- 1
    vline_type = "solid"
  }

  p <- ggplot()+
    geom_vline(xintercept = z0, color="grey80")+
    geom_density(aes(z0), color="black")+
    # geom_text(aes(x = (quantile(z0, 0.5) + z1)/2, y = 0.0003,
    #               label = sprintf("Monte-carlo, p=%.02f",p_value)),
    #           size = 3
    # ) +
    #geom_vline(xintercept = qnorm(c(0.05, 0.95),mean(z0),sd(z0)), color="blue")+
    geom_vline(xintercept = z1,color=vline_col, size = vline_size, linetype = vline_type)+
    #ggprism::theme_prism()+
    theme_linedraw()+
    theme(axis.text  = element_text(size = 10, face="bold"),
          axis.title  = element_text(size = 10, face="bold"),
          plot.title = element_text(size=12, face = "bold")) +
    theme(panel.grid.major = element_line(linetype = "blank"),
          panel.grid.minor = element_line(linetype = "blank"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()
    )+
    ggtitle(label = waiver(), subtitle = nonb_type)

  if (p_value < 0.001) {
    pwithsig <- p+
      geom_text(aes(x = (quantile(z0, 0.5) + z1)/2, y = 0.0003,
                    label = "p < 0.001, ***"),
                size = 3
      )
  } else if (p_value >= 0.05) {
    pwithsig <- p+
      geom_text(aes(x = (quantile(z0, 0.5) + z1)/2, y = 0.0003,
                    label = sprintf("NS")),
                size = 3
      )
  } else {
    pwithsig <- p+
      geom_text(aes(x = (quantile(z0, 0.5) + z1)/2, y = 0.0003,
                    label = "p < 0.05, *"),
                size = 3
      )
  }

  return(p)

}
PlotVlineBatch <- function(all_nonB, gtitle=NULL ,
                           metric = "total_overlapped_region", ...) {

  nonb_names <- unique(all_nonB$motif_type)


  dist_list <- map(.x = nonb_names, .f = PlotVline_clean,
                   all_nonB = all_nonB, metric = metric,
                   x_title= switch(
                     metric,
                     ho_rate = "Embeddability (average motifs per region)",
                     total_overlapped_region = "Oxidative regions (co-localized)",
                     total_overlapped_motif =  "Number of nonB (co-localized)"
                   ))
  patchwork::wrap_plots(dist_list, ...)+
    patchwork::plot_annotation(title = gtitle)
}
