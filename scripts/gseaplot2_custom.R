gseaplot2_custom<- function(
    x,
    geneSetID,
    title = "",
    color = "green",
    base_size = 11,
    rel_heights = c(1.5, .5, 1),
    subplots = 1:3,
    pvalue_table = FALSE,
    pvalue_table_columns = c("pvalue", "p.adjust"),
    pvalue_table_rownames = "Description",
) {
  ES_geom <- "line"
  
  geneList <- position <- NULL ## to satisfy codetool
  
  gsdata <- get_gsdata(x, geneSetID)
  
  p <- ggplot(gsdata, aes(x = .data$x)) +
    xlab(NULL) +
    theme_classic(base_size) +
    theme(
      panel.grid.major = element_line(colour = "grey92"),
      panel.grid.minor = element_line(colour = "grey92"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0))
  
  if (ES_geom == "line") {
    es_layer <- geom_line(
      aes(y = .data$runningScore, color = .data$Description),
      linewidth = 1
    )
  } else {
    es_layer <- geom_point(
      aes(y = .data$runningScore, color = .data$Description),
      size = 1,
      data = subset(gsdata, position == 1)
    )
  }
  
  p.res <- p +
    es_layer +
    theme(
      legend.position = "inside",
      legend.position.inside = c(.8, .8),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent")
    )
  
  p.res <- p.res +
    ylab("Running Enrichment Score") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(t = .2, r = .2, b = 0, l = .2, unit = "cm")
    )
  
  # Vectorized ymin/ymax assignment
  terms <- unique(gsdata$Description)
  term_indices <- match(gsdata$Description, terms) - 1
  idx <- which(gsdata$ymin != 0)
  gsdata[idx, "ymin"] <- term_indices[idx]
  gsdata[idx, "ymax"] <- term_indices[idx] + 1
  p2 <- ggplot(gsdata, aes(x = .data$x)) +
    geom_linerange(aes(
      ymin = .data$ymin,
      ymax = .data$ymax,
      color = .data$Description
    )) +
    xlab(NULL) +
    ylab(NULL) +
    theme_classic(base_size) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = -.1, b = 0, unit = "cm"),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line.x = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  
  if (length(geneSetID) == 1) {
    ## geneList <- gsdata$geneList
    ## j <- which.min(abs(geneList))
    ## v1 <- quantile(geneList[1:j], seq(0,1, length.out=6))[1:5]
    ## v2 <- quantile(geneList[j:length(geneList)], seq(0,1, length.out=6))[1:5]
    
    ## v <- sort(c(v1, v2))
    ## inv <- findInterval(geneList, v)
    
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) {
      inv <- inv + 1
    }
    
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(
      ymin = ymin,
      ymax = yy,
      xmin = xmin,
      xmax = xmax,
      col = col[unique(inv)]
    )
    p2 <- p2 +
      geom_rect(
        aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = .data$ymin,
          ymax = .data$ymax,
          fill = I(col)
        ),
        data = d,
        alpha = .9,
        inherit.aes = FALSE
      )
  }
  
  ## p2 <- p2 +
  ## geom_rect(aes(xmin=x-.5, xmax=x+.5, fill=geneList),
  ##           ymin=ymin, ymax = ymin + yy, alpha=.5) +
  ## theme(legend.position="none") +
  ## scale_fill_gradientn(colors=color_palette(c("blue", "red")))
  
  df2 <- p$data #data.frame(x = which(p$data$position == 1))
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p +
    geom_segment(
      data = df2,
      aes(x = .data$x, xend = .data$x, y = .data$y, yend = 0),
      color = "grey"
    )
  p.pos <- p.pos +
    ylab("Ranked List Metric") +
    xlab("Rank in Ordered Dataset") +
    theme(
      plot.margin = margin(t = -.1, r = .2, b = .2, l = .2, unit = "cm")
    )
  
  if (!is.null(title) && !is.na(title) && title != "") {
    p.res <- p.res + ggtitle(title)
  }
  
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    } else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  
  if (pvalue_table) {
    pd <- x[geneSetID, pvalue_table_columns]
    # pd <- pd[order(pd[,1], decreasing=FALSE),]
    if (is.null(pvalue_table_rownames)) {
      rows <- NULL
    } else {
      # rownames(pd) <- pd$Description
      if (length(pvalue_table_rownames) != 1) {
        stop(
          "the length of `pvalue_table_rownames` should be equal to 1"
        )
      }
      
      rows <- x[geneSetID, pvalue_table_rownames]
    }
    
    # pd <- round(pd, 4)
    for (i in seq_len(ncol(pd))) {
      pd[, i] <- format(pd[, i], digits = 4)
    }
    tp <- tableGrob2(d = pd, p = p.res, rows = rows)
    
    p.res <- p.res +
      theme(legend.position = "none") +
      annotation_custom(
        tp,
        xmin = quantile(p.res$data$x, .5),
        xmax = quantile(p.res$data$x, .95),
        ymin = quantile(p.res$data$runningScore, .75),
        ymax = quantile(p.res$data$runningScore, .9)
      )
  }
  
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(
      axis.line.x = element_line(),
      axis.ticks.x = element_line(),
      axis.text.x = element_text()
    )
  
  if (length(subplots) == 1) {
    return(
      plotlist[[1]] +
        theme(
          plot.margin = margin(
            t = .2,
            r = .2,
            b = .2,
            l = .2,
            unit = "cm"
          )
        )
    )
  }
  
  if (length(rel_heights) > length(subplots)) {
    rel_heights <- rel_heights[subplots]
  }
  
  # aplot::plot_list(gglist = plotlist, ncol=1, heights=rel_heights)
  aplot::gglist(gglist = plotlist, ncol = 1, heights = rel_heights)
}


gsInfo_custom <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID)) {
    geneSetID <- object@result[geneSetID, "ID"]
  }
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  } else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}