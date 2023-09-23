##' extract gsea result of selected geneSet
##'
##'
##' @title gsInfo
##' @param object gseaResult object
##' @param geneSetID gene set ID
##' @return data.frame
##' @author Guangchuang Yu
## @export
gseaScores <- utils::getFromNamespace("gseaScores", "DOSE")
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

## Modify gseaplot2 function ###

gseaplot2_mod_tt <- function(x, geneSetID, title = "", color="green", color_hmp = c("#39518F","#7D2125"), base_size = 11, fontface = "plain",
                      rel_heights=c(1.5, .5, 1), subplots = 1:3, nesDigit = 2, pDigit = 2, pvalSize = 4,pCol = "grey30", pHjust = 1,
                      pvalue_table = FALSE, addPval = TRUE, pvalX = 0.9, pvalY = 0.9,  ES_geom="line") {
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  
  geneList <- position <- NULL ## to satisfy codetool
  
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  } else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(base_size) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand=c(0,0))
  
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description),
                          size=1)
  } else {
    es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description),
                           size=1, data = subset(gsdata, position == 1))
  }
  
  p.res <- p + es_layer +
    theme(legend.position = c(.8, .8), legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"))
  
  p.res <- p.res + ylab("Running Enrichment Score") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          #text=element_text(family="Arial"),
          plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
  
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
    theme(legend.position = "none",
          #text=element_text(family="Arial"),
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out=9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) inv <- inv + 1
    
    col <- grDevices::colorRampPalette(c(color_hmp[1], "white", color_hmp[2]))(10)  #c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy,
                    xmin = xmin,
                    xmax = xmax,
                    col = col[unique(inv)])
    p2 <- p2 + geom_rect(
      aes_(xmin=~xmin,
           xmax=~xmax,
           ymin=~ymin,
           ymax=~ymax,
           fill=~I(col)),
      data=d,
      alpha=.9,
      inherit.aes=FALSE)
  }
  
  df2 <- p$data #data.frame(x = which(p$data$position == 1))
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                            color="grey")
  p.pos <- p.pos + ylab("Ranked List Metric") +
    xlab("Rank in Ordered Dataset") +
    theme(#text=element_text(family="Arial"), 
          plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
  
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)
  
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values=color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    } else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    # pd <- pd[order(pd[,1], decreasing=FALSE),]
    rownames(pd) <- pd$Description
    
    pd <- pd[,-1]
    # pd <- round(pd, 4)
    for (i in seq_len(ncol(pd))) {
      pd[, i] <- format(pd[, i], digits = 4)
    }
    tp <- tableGrob2(pd, p.res)
    
    p.res <- p.res + theme(legend.position = "none") +
      annotation_custom(tp,
                        xmin = quantile(p.res$data$x, .5),
                        xmax = quantile(p.res$data$x, .95),
                        ymin = quantile(p.res$data$runningScore, .75),
                        ymax = quantile(p.res$data$runningScore, .9))
  }
  
  # to dataframe
  data_ga <- data.frame(x) %>%
    dplyr::filter(ID %in% geneSetID)
  data_ga <- data_ga[unique(gsdata$Description),]
  
  # add NES Pvalue
  if (addPval == TRUE) {
    pLabel <- paste0(
      "NES = ",
      round(data_ga$NES, digits = nesDigit),
      "\n",
      "Pvalue = ",
      # round(data_ga$pvalue, digits = pDigit),
      ifelse(data_ga$pvalue < 0.001,"< 0.001",round(data_ga$pvalue, digits = pDigit)),
      "\n",
      "Adjusted Pvalue = ",
      # round(data_ga$p.adjust, digits = pDigit),
      ifelse(data_ga$p.adjust < 0.001,"< 0.001",round(data_ga$p.adjust, digits = pDigit)),
      "\n",
      sep = " "
    )
    
    px <- pvalX * nrow(gsdata[which(gsdata$Description == geneSetID[1]),])
    py <-
      pvalY * sum(abs(range(gsdata$runningScore))) + min(gsdata$runningScore)
    
    # add pvlaue label
    
    if(length(geneSetID) == 1){
      p.res <-
        p.res +
        ggplot2::annotate(geom = "text",
                          x = px,
                          y = py,
                          label = pLabel,
                          size = pvalSize,
                          color = pCol,
                          fontface = fontface,
                          #family = "Arial",
                          hjust = pHjust)
    }else{
      mytable <- tibble::tibble(x = px, y = py,
                                table = list(tibble::tibble('NES' = round(data_ga$NES, digits = nesDigit),
                                                            'Pvalue' = ifelse(data_ga$pvalue < 0.001,"< 0.001",round(data_ga$pvalue, digits = pDigit)),
                                                            'Adjusted Pvalue' = ifelse(data_ga$p.adjust < 0.001,"< 0.001",round(data_ga$p.adjust, digits = pDigit)))))
      
      
      p.res <- p.res +
        ggpp::geom_table(data = mytable, ggplot2::aes(px, py, label = table))
    }
    
  } else {
    p.res <- p.res
  }
  
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x=element_line(),
          axis.text.x = element_text())
  
  
  if (length(subplots) == 1)
    return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2,
                                                    l=.2, unit="cm")))
  
  
  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]
  
  # aplot::plot_list(gglist = plotlist, ncol=1, heights=rel_heights)
  pfinal = aplot::gglist(gglist = plotlist, ncol=1, heights=rel_heights) 
  
  return(pfinal)
}
