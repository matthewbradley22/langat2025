idents.use = 'Astrocytes'
color.use = c("grey10", "#F8766D", "#00BFC4")
comparison = c(1,2)
signaling = NULL
signaling.label = NULL
top.label = 1
signaling.exclude = NULL
xlims = NULL
ylims = NULL
slot.name = "netP"
dot.size = 2.5
point.shape = c(21, 22, 24, 23)
label.size = 3
dot.alpha = 0.6
x.measure = "outdeg"
y.measure = "indeg"
xlabel = "Differential outgoing interaction strength"
ylabel = "Differential incoming interaction strength"
title = NULL
font.size = 10
font.size.title = 10
do.label = T
show.legend = T
show.axes = T

#Taken from netAnalysis_signalingChanges_scatter function on cellchat github
netAnalysis_signalingChanges_custom <- function(object, return_obj = TRUE){
  if (is.list(object)) {
    object <- mergeCellChat(object, add.names = names(object))
  }
  if (is.list(object@net[[1]])) {
    dataset.name <- names(object@net)
    message(paste0("Visualizing differential outgoing and incoming signaling changes from ", dataset.name[comparison[1]], " to ", dataset.name[comparison[2]]))
    title <- paste0("Signaling changes of ", idents.use, " (", dataset.name[comparison[1]], " vs. ", dataset.name[comparison[2]], ")")
    
    cell.levels <- levels(object@idents$joint)
    if (is.null(xlabel) | is.null(ylabel)) {
      xlabel = "Differential outgoing interaction strength"
      ylabel = "Differential incoming interaction strength"
    }
    
  } else {
    message("Visualizing outgoing and incoming signaling on a single object \n")
    title <- paste0("Signaling patterns of ", idents.use)
    if (length(slot(object, slot.name)$centr) == 0) {
      stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
    }
    cell.levels <- levels(object@idents)
  }
  if (!(idents.use %in% cell.levels)) {
    stop("Please check the input cell group names!")
  }
  if (is.null(signaling)) {
    signaling <- union(object@netP[[comparison[1]]]$pathways, object@netP[[comparison[2]]]$pathways)
  }
  if (!is.null(signaling.exclude)) {
    signaling <- setdiff(signaling, signaling.exclude)
  }
  mat.all.merged <- list()
  for (ii in 1:length(comparison)) {
    if (length(slot(object, slot.name)[[comparison[ii]]]$centr) == 0) {
      stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores for each dataset seperately! ")
    }
    if (sum(c(x.measure, y.measure) %in% names(slot(object, slot.name)[[comparison[ii]]]$centr[[1]])) !=2) {
      stop(paste0("`x.measure, y.measure` should be one of ", paste(names(slot(object, slot.name)[[comparison[ii]]]$centr[[1]]),collapse=", "), '\n', "`outdeg_unweighted` is only supported for version >= 1.1.2"))
    }
    centr <- slot(object, slot.name)[[comparison[ii]]]$centr
    outgoing <- matrix(0, nrow = length(cell.levels), ncol = length(centr))
    incoming <- matrix(0, nrow = length(cell.levels), ncol = length(centr))
    dimnames(outgoing) <- list(cell.levels, names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
      outgoing[,i] <- centr[[i]][[x.measure]]
      incoming[,i] <- centr[[i]][[y.measure]]
    }
    mat.out <- t(outgoing)
    mat.in <- t(incoming)
    
    mat.all <- array(0, dim = c(length(signaling),ncol(mat.out),2))
    mat.t <-list(mat.out, mat.in)
    for (i in 1:length(comparison)) {
      mat = mat.t[[i]]
      mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
      mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
      idx <- match(rownames(mat1), signaling)
      mat[idx[!is.na(idx)], ] <- mat1
      dimnames(mat) <- list(signaling, colnames(mat1))
      mat.all[,,i] = mat
    }
    dimnames(mat.all) <- list(dimnames(mat)[[1]], dimnames(mat)[[2]], c("outgoing", "incoming"))
    mat.all.merged[[ii]] <- mat.all
  }
  mat.all.merged.use <- list(mat.all.merged[[1]][,idents.use,], mat.all.merged[[2]][,idents.use,])
  idx.specific <- mat.all.merged.use[[1]] * mat.all.merged.use[[2]]
  mat.sum <- mat.all.merged.use[[2]] +  mat.all.merged.use[[1]]
  out.specific.signaling <- rownames(idx.specific)[(mat.sum[,1] != 0) & (idx.specific[,1] == 0)]
  in.specific.signaling <- rownames(idx.specific)[(mat.sum[,2] != 0) & (idx.specific[,2] == 0)]
  
  mat.diff <- mat.all.merged.use[[2]] -  mat.all.merged.use[[1]]
  idx <- rowSums(mat.diff) != 0
  mat.diff <- mat.diff[idx, ]
  out.specific.signaling <- rownames(mat.diff) %in% out.specific.signaling
  in.specific.signaling <- rownames(mat.diff) %in% in.specific.signaling
  out.in.specific.signaling <- as.logical(out.specific.signaling * in.specific.signaling)
  specificity.out.in <- matrix(0, nrow = nrow(mat.diff), ncol = 1)
  specificity.out.in[out.in.specific.signaling] <- 2 # both outgoing and incoming specific to one condition
  specificity.out.in[setdiff(which(out.specific.signaling), which(out.in.specific.signaling))] <- 1 # only outgoing specific to one condition
  specificity.out.in[setdiff(which(in.specific.signaling), which(out.in.specific.signaling))] <- -1 # only incoming specific to one condition
  
  
  df <- as.data.frame(mat.diff)
  df$specificity.out.in <- specificity.out.in
  df$specificity = 0
  df$specificity[(specificity.out.in != 0) & (rowSums(mat.diff >= 0) ==2)] = 1 # specific to dataset 2
  df$specificity[(specificity.out.in != 0) & (rowSums(mat.diff <= 0) ==2)] = -1  # specific to dataset 1
  
  # change number to char
  out.in.category <- c("Shared", "Incoming specific", "Outgoing specific", "Incoming & Outgoing specific")
  specificity.category <- c("Shared", paste0(dataset.name[comparison[1]]," specific"), paste0(dataset.name[comparison[2]]," specific"))
  df$specificity.out.in <- plyr::mapvalues(df$specificity.out.in, from = c(0,-1,1,2),to = out.in.category)
  df$specificity.out.in <- factor(df$specificity.out.in, levels = out.in.category)
  df$specificity <- plyr::mapvalues(df$specificity, from = c(0,-1,1),to = specificity.category)
  df$specificity <- factor(df$specificity, levels = specificity.category)
  
  point.shape.use <- point.shape[out.in.category %in% unique(df$specificity.out.in)]
  df$specificity.out.in = droplevels(df$specificity.out.in, exclude = setdiff(out.in.category,unique(df$specificity.out.in)))
  
  color.use <- color.use[specificity.category %in% unique(df$specificity)]
  df$specificity = droplevels(df$specificity, exclude = setdiff(specificity.category,unique(df$specificity)))
  
  df$labels <- rownames(df)
  
  if(return_obj == TRUE){
    return(df)
  }
  gg <- ggplot(data = df, aes(outgoing, incoming)) +
    geom_point(aes(colour = specificity, fill = specificity, shape = specificity.out.in), size = dot.size)
  gg <- gg + theme_linedraw() +theme(panel.grid = element_blank()) +
    geom_hline(yintercept=0,linetype="dashed", color = "grey50", size = 0.25) + geom_vline(xintercept=0, linetype="dashed", color = "grey50",size = 0.25) +
    theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15, "in"))+
    # guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) + theme(plot.title = element_text(size= font.size.title, hjust = 0.5, face="plain"))+
    # theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE) + guides(fill="none")
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
  gg <- gg + scale_shape_manual(values = point.shape.use)
  gg <- gg + theme(legend.title = element_blank())
  if (!is.null(xlims)) {
    gg <- gg + xlim(xlims)
  }
  if (!is.null(ylims)) {
    gg <- gg + ylim(ylims)
  }
  
  if (do.label) {
    if (is.null(signaling.label)) {
      thresh <- stats::quantile(abs(as.matrix(df[,1:2])), probs = 1-top.label)
      idx = abs(df[,1]) > thresh | abs(df[,2]) > thresh
      data.label <- df[idx,]
    } else {
      data.label <- df[rownames(df) %in% signaling.label, ]
    }
    
    gg <- gg + ggrepel::geom_text_repel(data = data.label, mapping = aes(label = labels, colour = specificity), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5)
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  
  if (!show.axes) {
    gg <- gg + theme_void()
  }
  
  gg
}

#Data for function created in sc_cellchat.R script
df_wt_3 <- netAnalysis_signalingChanges_custom(cellchat_wt_3_merged)
df_wt_4 <- netAnalysis_signalingChanges_custom(cellchat_wt_4_merged)

df_ips_3 <- netAnalysis_signalingChanges_custom(cellchat_ips_3_merged)
df_ips_4 <- netAnalysis_signalingChanges_custom(cellchat_ips_4_merged)
df_ips_5 <- netAnalysis_signalingChanges_custom(cellchat_ips_5_merged)

#Which pathways have large differences between datasets
df_wt_3 %>%  dplyr::arrange(desc(incoming))
df_ips_3 %>%  dplyr::arrange(desc(incoming))

#Plot specific pathways of interest
netVisual_bubble(cellchat_wt_3_merged, targets.use = 'Astrocytes',  comparison = c(1, 2), angle.x = 45, signaling = 'Glutamate')
netVisual_bubble(cellchat_ips_3_merged, targets.use = 'Astrocytes',  comparison = c(1, 2), angle.x = 45, signaling = 'SEMA4')

#Plot incoming pathways that increase in ips

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/wt_3_vs_mock_up_cellchat.pdf', height = 6, width = 6)
df_wt_3 %>%  dplyr::arrange(desc(incoming)) %>% 
  head(n = 15) %>% 
  ggplot(aes(x = incoming, y = reorder(labels, incoming), fill = specificity))+
  geom_bar(stat = 'identity', color = 'black')+
  theme_classic()+
  ylab('')+
  xlab('Incoming interaction strength change')+
  theme(text = element_text(size = 18))+
  xlim(c(0, 2.9))+
  scale_fill_manual(values = c("#CFB3E2", "#B4DBDB"))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/wt_4_vs_mock_up_cellchat.pdf', height = 6, width = 6)
df_wt_4%>%  dplyr::arrange(desc(incoming)) %>% 
  head(n = 15) %>% 
  ggplot(aes(x = incoming, y = reorder(labels, incoming), fill = specificity))+
  geom_bar(stat = 'identity', color = 'black')+
  theme_classic()+
  ylab('')+
  xlab('Incoming interaction strength change')+
  theme(text = element_text(size = 18))+
  xlim(c(0, 2.9))+
  scale_fill_manual(values = c("#CFB3E2", "#B4DBDB"))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/ips_3_vs_mock_up_cellchat.pdf', height = 6, width = 6)
df_ips_3 %>%  dplyr::arrange(desc(incoming)) %>% 
  head(n = 15) %>% 
  ggplot(aes(x = incoming, y = reorder(labels, incoming), fill = specificity))+
  geom_bar(stat = 'identity', color = 'black')+
  theme_classic()+
  ylab('')+
  xlab('Incoming interaction strength change')+
  theme(text = element_text(size = 18))+
  xlim(c(0, 2.9))+
  scale_fill_manual(values = c("#CFB3E2", "#B4DBDB"))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/ips_4_vs_mock_up_cellchat.pdf', height = 6, width = 6)
df_ips_4 %>%  dplyr::arrange(desc(incoming)) %>% 
  head(n = 15) %>% 
  ggplot(aes(x = incoming, y = reorder(labels, incoming), fill = specificity))+
  geom_bar(stat = 'identity', color = 'black')+
  theme_classic()+
  ylab('')+
  xlab('Incoming interaction strength change')+
  theme(text = element_text(size = 18))+
  xlim(c(0, 2.9))+
  scale_fill_manual(values = c("#CFB3E2", "#B4DBDB"))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/ips_5_vs_mock_up_cellchat.pdf', height = 6, width = 6)
df_ips_5 %>%  dplyr::arrange(desc(incoming)) %>% 
  head(n = 15) %>% 
  ggplot(aes(x = incoming, y = reorder(labels, incoming), fill = specificity))+
  geom_bar(stat = 'identity', color = 'black')+
  theme_classic()+
  ylab('')+
  xlab('Incoming interaction strength change')+
  theme(text = element_text(size = 18))+
  xlim(c(0, 2.9))+
  scale_fill_manual(values = c("#CFB3E2", "#B4DBDB"))
dev.off()
