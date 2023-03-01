##########################
#' Plot covariance matrices as a network
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @details This function is meant for internal use only.
#' @export
#' 
get_cov_plot <- function(indexes, CORRS, DIAG, A){
  
  std_lim <- range(sapply(indexes, function(ii) sqrt(range(DIAG[[ii]]))))
  
  plots <- list()
  count <- 1
  for(kkk in indexes){
    CO <- CORRS[[kkk]] #- CORRS[[which.min(aver_corr)]]
    STD <- sqrt(DIAG[[kkk]])
    
    nel <- nrow(A)
    xp <- cos(2*pi*(1:nel/nel))
    yp <- sin(2*pi*(1:nel/nel))
    
    if(nrow(A) == 14){
      xp <- c(0, 0, 0.2, 0.2, -0.1, -0.2, -0.2, 0, 0.2, 0.35, 0.2, 0.35, 0, -0.2)
      yp <- c(0.9, 0.7, 0.6, 0.3, 0.5, 0.1, -0.1, 0, 0.05, -0.1, -0.3, -0.45, -0.46, -0.5)
    } else {
      xp <- c(0.2, 0.21, 0.25, 0.3, 0.17)
      yp <- c(0.3, 0.13, 0, -0.1, -0.2)
    }
    
    dat <- data.frame("x" = xp, "y" = yp, "std" = STD,
                      "lab" = if(nrow(A)==14){ 1:length(STD) } else {c("Sco", "Nor", "Mid", "Lon", "Sou")})
    dat_cor <- list()
    kk <- 1
    for(ir in 1:(nel-1))
      for(ic in (ir+1):nel){
        dat_cor[[kk]] <- c(xp[ir],yp[ir],xp[ic],yp[ic],CO[ir, ic])
        kk <- kk + 1
      }
    
    dat_cor <- as.data.frame(do.call("rbind", dat_cor))
    colnames(dat_cor) <- c("x", "y", "xend", "yend", "cor")
    dat_cor_09 <- dat_cor[dat_cor$cor > 0.9, ]
    #dat_cor <- rbind(dat_cor, dat_cor_0.9)
    
    if(nrow(A)==14){
      stdleg <- 1
      corleg <- 2
    } else {
      stdleg <- 2
      corleg <- Inf
    }
    library(ggplot2)
    pl <- ggplot(data = dat, mapping = aes(x = x, y = y)) +
      geom_segment(data = dat_cor, mapping = aes(x = x, y = y, xend = xend, yend = yend,
                                                 colour = cor, alpha = abs(cor)),
                   inherit.aes = FALSE, size = 2) +
      geom_point(mapping = aes(x = x, y = y, fill = std),
                 size = if(nrow(A)==14){20}else{20}, shape = 21) +
      geom_text(aes(label = lab), size = if(nrow(A)==14){10}else{8}, colour = "white") +
      scale_fill_viridis_c(option = "magma", end = 0.8, name = "Std. dev.", limits = std_lim,
                           guide=if(count==stdleg){ guide_colorbar() } else { "none" }) +
      scale_colour_gradient(low="darkblue", high="red", limits = c(0, 1), name = "Correlation",
                            guide=if(count==corleg){ guide_colorbar() } else { "none" }) +
      #scale_colour_viridis_c(option = "magma", end = 0.9, name = "Corr") +
      scale_alpha_continuous(limits = c(0, 1),guide="none") +
      #scale_alpha_continuous(range = c(0, 1))  +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            panel.border = element_blank(),
            legend.position = "bottom")
    if(nrow(A)!=14){
      pl <- pl + xlim(0.15, 0.31) + ylim(-0.3, 0.33)
   #   xp <- c(0.2, 0.21, 0.25, 0.3, 0.17)
    #  yp <- c(0.3, 0.13, 0, -0.1, -0.2)
    }
    plots[[count]] <- pl
    count <- count + 1
    
  }
  
  return(plots)
  
}
