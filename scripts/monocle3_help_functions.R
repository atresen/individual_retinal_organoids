## functions compatible with monocle3 ##

plot_umap_gene_expr <- function(cds, marker, cell_size = 1, stroke_size = 0.0){
                           
  cds_subset <- cds[rowData(cds)$gene_short_name == marker,]
  marker_expr <- Matrix::t(assay(cds_subset))/colData(cds_subset)$Size_Factor
  
  meta_df <- as.data.frame(colData(cds_subset))
  meta_df$umap1 <- reducedDim(x = cds, type = "UMAP")[,1]
  meta_df$umap2 <- reducedDim(x = cds, type = "UMAP")[,2]
  
  meta_df$marker <- as.numeric(marker_expr[,1])
  meta_df$marker = with(meta_df, ifelse(marker >= 0.01, marker, NA))
  
  ggplot(meta_df, aes(x = umap1, y = umap2, color = log10(marker + 0.1))) +
    geom_point(size = cell_size, stroke = stroke_size, aes(alpha = ifelse(!is.na(marker), "1", "2"))) + 
    scale_color_viridis(option = "viridis",
          name = "log10(",marker," + 0.1)", na.value = "grey80") +
    scale_alpha_manual(values = c(1, .3)) + guides(alpha = FALSE) +
    monocle3:::monocle_theme_opts() +
    theme_void()
}


# same function as above but modified so that the legend is plotted too
plot_umap_gene_expr_legend = function (cds, marker, cell.size = 0.5) {
    
    cds_subset <- cds[rowData(cds)$gene_short_name == marker, ]
    marker_expr <- Matrix::t(assay(cds_subset))/colData(cds_subset)$Size_Factor

    meta_df <- as.data.frame(colData(cds_subset))
  	meta_df$umap1 <- reducedDim(x = cds, type = "UMAP")[,1]
    meta_df$umap2 <- reducedDim(x = cds, type = "UMAP")[,2]

    meta_df$marker <- as.numeric(marker_expr[, 1])
    meta_df$marker = with(meta_df, ifelse(marker >= 0.01, marker, NA))

    color_str <- paste0("log10(",marker," + 0.1)")
    ggplot(meta_df, aes(x = umap1, y = umap2, color = log10(marker +
        0.1))) + geom_point(size = cell.size, stroke = 0.1, aes(alpha = ifelse(!is.na(marker),
        "1", "2"))) + scale_color_viridis(option = "viridis",
        name = color_str, na.value = "grey80") + 
		scale_alpha_manual(values = c(1, 0.5)) +
        guides(alpha = FALSE) +
        theme_void() +
        theme(legend.position = "right")
}


# plot overlapping expression of 2 genes - this is the umap version
plot.two.genes = function(cds, gene1, gene2, thresh = 1) {
    gene.id1 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene1, "id"])
    gene.id2 = as.character(rowData(cds)[rowData(cds)$gene_short_name == gene2, "id"]) 
    
    colData(cds)$umap_1 <- reducedDim(x = cds, type = "UMAP")[,1]
	 colData(cds)$umap_2 <- reducedDim(x = cds, type = "UMAP")[,2]

    if(sum(assay(cds)[gene.id1,] >= thresh & assay(cds)[gene.id2,] >= thresh) > 0) {
        colData(cds)$tmp = ifelse(
                assay(cds)[gene.id1,] >= thresh & assay(cds)[gene.id2,] >= thresh, "both",
                ifelse(assay(cds)[gene.id1,] >= thresh, gene1, 
                      ifelse(assay(cds)[gene.id2,] >= thresh, gene2, "neither")))
        colData(cds)$tmp = factor(colData(cds)$tmp, levels = c(gene1, gene2, "both", "neither"))
        plot = ggplot(as.data.frame(colData(cds)), aes(x = umap_1, y = umap_2, color = tmp)) +
            geom_point(size = 0.2, aes(
                alpha = ifelse(assay(cds)[gene.id1,] > 0 | 
                               assay(cds)[gene.id2,] > 0, "1", "2"))) +
            scale_color_manual(values = c("royalblue3", "firebrick3", "mediumorchid3" , "grey80")) +
            scale_alpha_manual(values = c(1, .5)) + guides(alpha = FALSE) +
            guides(colour = guide_legend(override.aes = list(size = 2))) +
            monocle3:::monocle_theme_opts() +
            theme_void()
    } else {
        colData(cds)$tmp = 
                    ifelse(assay(cds)[gene.id1,] >= thresh, gene1, 
                    ifelse(assay(cds)[gene.id2,] >= thresh, gene2, "neither"))
        colData(cds)$tmp = factor(colData(cds)$tmp, levels = c(gene1, gene2, "neither"))
        plot = ggplot(as.data.frame(colData(cds)), aes(x = umap_1, y = umap_2, color = tmp)) +
            geom_point(size = 0.2, aes(
                alpha = ifelse(assay(cds)[gene.id1,] > 0 | 
                               assay(cds)[gene.id2,] > 0, "1", "2"))) +
            scale_color_manual(values = c("royalblue3", "firebrick3", "grey80")) +
            scale_alpha_manual(values = c(1, .5)) + guides(alpha = FALSE) +
            guides(colour = guide_legend(override.aes = list(size = 2))) +
            monocle3:::monocle_theme_opts() +
            theme_void()
    }
    
    colData(cds)$tmp = NULL
    return(plot)
}



# generate an array of expression plots
save_umap_expr_array <- function(cds, marker_genes, file = "tmp_marker_array.png", n.row = 1,
                                  plot.width = 4, plot.height = 4, cell.size = 1){
  
  plot.list <- list()
  
  for(marker in marker_genes){
    cds_subset <- cds[rowData(cds)$gene_short_name == marker, ]
    marker_expr <- Matrix::t(assay(cds_subset))/colData(cds_subset)$Size_Factor
    
    meta_df <- as.data.frame(colData(cds_subset))
    meta_df$umap1 = reducedDim(x = cds, type = "UMAP")[,1]
    meta_df$umap2 = reducedDim(x = cds, type = "UMAP")[,2]
    
    meta_df$marker <- as.numeric(marker_expr[, 1])
    meta_df$marker = with(meta_df, ifelse(marker >= 0.01, marker, NA))
    
    color_str <- paste0("log10(",marker," + 0.1)")
    
    plot.list[[marker]] <- ggplot(meta_df, aes(x = umap1, y = umap2, color = log10(marker + 0.1))) +
      geom_point(size = cell.size, stroke = 0.0, aes(alpha = ifelse(!is.na(marker), "1", "2"))) + 
      ggtitle(marker) +
      scale_color_viridis(option = "viridis", direction = 1,
                          name = "log10(",marker," + 0.1)", na.value = "grey75") +
      theme_void() + scale_alpha_manual(values = c(1, .3)) + guides(alpha = FALSE) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "right", legend.title = element_text(size = 6), 
            legend.text = element_text(size = 6), legend.key.size = unit(.01, "npc"), 
            legend.key.height = unit(.02, "npc"), text = element_text(size = 6))
    
  }
  
  png(file, units = "in", width = plot.width, height = plot.height, res = 600)
  do.call("grid.arrange", c(plot.list, nrow = n.row))
  dev.off()
}


# calculate signature scores for groups of genes (from Jose)
calculate_signature_scores = function (cds, gene_list) {
    vst = as.data.frame(as.matrix(monocle3::normalized_counts(cds, norm_method = "log")))
    vst = vst[gene_list$id, ]
    vst = merge(vst, gene_list[, names(gene_list) %in% c("id", 
        "signature")], by.x = "row.names", by.y = "id")
    vst <- melt(vst, id.vars = c("Row.names", "signature"), variable.name = "cell", 
        value.name = "vst")
    scores = vst %>% group_by(cell, signature) %>% summarize(score = mean(vst))
    final_scores = recast(scores, signature ~ cell)
    final_scores = t(final_scores)
    colnames(final_scores) = final_scores[1, ]
    final_scores = as.data.frame(final_scores[-1, ])
    for (column in colnames(final_scores)) {
        final_scores[, column] <- as.numeric(as.character(final_scores[, 
            column]))
    }
    final_scores = scale(final_scores)
    return(final_scores)
}

# calculate scores (updated version)
estimate_gene_scores = function(cds, gene_list1, gene_list2){
  
  cds1 = cds[fData(cds)$gene_short_name %in% gene_list1,]
  aggregate_expression_1 = assay(cds1)
  aggregate_expression_1 = t(t(aggregate_expression_1) / colData(cds1)$Size_Factor)
  aggregate_expression_1 = Matrix::colSums(aggregate_expression_1)
  
  cds_g2m = cds[fData(cds)$gene_short_name %in% gene_list2,]
  aggregate_g2m_expression = assay(cds_g2m)
  aggregate_g2m_expression = t(t(aggregate_g2m_expression) / pData(cds_g2m)$Size_Factor)
  aggregate_g2m_expression = Matrix::colSums(aggregate_g2m_expression)

  pData(cds)$g1s_score = log(aggregate_g1s_expression+1)
  pData(cds)$g2m_score = log(aggregate_g2m_expression+1)
  
  pData(cds)$proliferation_index = log(aggregate_g1s_expression + aggregate_g2m_expression + 1)
  
  return(cds)

  }

# return a list of gene ids from short names
get.gene.ids = function (cds, names) {
  return(as.character(rowData(cds)[rowData(cds)$gene_short_name %in% names, "id"]))
}

### generate 3D expression plots with plotly (from JP)
plot.expr.3d = function(cds, gene_short_name, file = NULL,
                        coexpr = F, size = 2,
                        color.pal = c("grey85", "red", "orange", "yellow"),
                        downsample = NULL) {
    
    if (coexpr & length(nchar(gene_short_name)) == 1) {
        message("Error. If coexpr = T, must give vector of >1 gene name as input.")
        return(NULL)
    }
    
    if (!is.null(downsample)) {
        set.seed(42)
        ids = sample(1:ncol(cds), downsample, replace = F)
    } else {
        ids = 1:ncol(cds)
    }
    
    if (is.null(file)) {
        this.dir = system("pwd", intern = T)
        
        if (coexpr) {
            file = paste(this.dir, "/plots/cil.", paste(gene_short_name, collapse = "."), ".coexpression.html", sep = "")
        } else {
            file = paste(this.dir, "/plots/cil.", gene_short_name, ".html", sep = "")
        }
    }
    
    if (coexpr) {
        tmp.df = rownames(colData(cds)) %>% as.data.frame() %>% select("cell" = ".")
        tmp.df$Cluster = clusters(cds)
        tmp.df$umap_1 = reducedDim(x = cds, type = "UMAP")[,1]
        tmp.df$umap_2 = reducedDim(x = cds, type = "UMAP")[,2]
        tmp.df$umap_3 = reducedDim(x = cds, type = "UMAP")[,3]
        gene.ids = sapply(gene_short_name, function(x) get.gene.id(cds, x))
        gene.sanitized.name = sub("-", ".", gene_short_name)

        for (i in 1:length(gene_short_name))
            tmp.df[, gene.sanitized.name[i]] = assay(cds)[gene.ids[i],] / colData(cds)$Size_Factor

        tmp.df$tmp.expr = log2(apply(tmp.df[, 6:(6+length(gene_short_name)-1)], 1, min)+1)
                          
        plot = plot_ly(
            tmp.df[ids,],
            x = ~ umap_1,
            y = ~ umap_2,
            z = ~ umap_3,
            marker = list(
                size = size,
                color = ~ tmp.expr,
                colorscale = color.pal,
                showscale = TRUE),
            text = ~ paste("Cluster: ", Cluster)) %>%
            add_markers()
                          
        htmlwidgets::saveWidget(plot, file)
    } else {
        gene.id = get.gene.id(cds, gene_short_name)
        colData(cds)$tmp.expr = log2(cds@assayData$exprs[gene.id,] / pData(cds)$Size_Factor + 1)

        plot = plot_ly(
            colData(cds)[ids,] %>% as.data.frame(),
            x = ~ umap_1,
            y = ~ umap_2,
            z = ~ umap_3,
            marker = list(
                size = size,
                color = ~ tmp.expr,
                colorscale = color.pal,
                showscale = TRUE),
            text = ~ paste("Cluster: ", Cluster, "<br>log2 norm expr: ", tmp.expr)) %>%
            add_markers()
        
        pData(cds)$tmp.expr = NULL
        
        htmlwidgets::saveWidget(plot, file)
    }
}

##### plot gene expression in 3D (from JP), this is for when there are no ensembl ids in the fData - just gene_short_names ######
plot.expr.3d.noid = function(cds, gene_short_name, file = NULL,
                        coexpr = F, size = 2,
                        color.pal = c("grey85", "red", "orange", "yellow"),
                        downsample = NULL) {
  
  if (coexpr & length(nchar(gene_short_name)) == 1) {
    message("Error. If coexpr = T, must give vector of >1 gene name as input.")
    return(NULL)
  }
  
  if (!is.null(downsample)) {
    set.seed(42)
    ids = sample(1:ncol(cds), downsample, replace = F)
  } else {
    ids = 1:ncol(cds)
  }
  
  if (is.null(file)) {
    this.dir = system("pwd", intern = T)
    
    if (coexpr) {
      file = paste(this.dir, "/plots/cil.", paste(gene_short_name, collapse = "."), ".coexpression.html", sep = "")
    } else {
      file = paste(this.dir, "/plots/cil.", gene_short_name, ".html", sep = "")
    }
  }
  
  if (coexpr) {
    tmp.df = rownames(colData(cds)) %>% as.data.frame() %>% select("cell" = ".")
        tmp.df$Cluster = clusters(cds)
        tmp.df$umap_1 = reducedDim(x = cds, type = "UMAP")[,1]
        tmp.df$umap_2 = reducedDim(x = cds, type = "UMAP")[,2]
        tmp.df$umap_3 = reducedDim(x = cds, type = "UMAP")[,3]
        gene.ids = sapply(gene_short_name, function(x) get.gene.id(cds, x))
    gene.sanitized.name = sub("-", ".", gene_short_name)
    
    for (i in 1:length(gene_short_name))
      tmp.df[, gene.sanitized.name[i]] = cds@assayData$exprs[gene_short_name[i],] / pData(cds)$Size_Factor
    
    tmp.df$tmp.expr = log2(apply(tmp.df[, 6:(6+length(gene_short_name)-1)], 1, min)+1)
    
    ax <- list(
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
    )

    plot = plot_ly(
      tmp.df[ids,],
      x = ~ umap_1,
      y = ~ umap_2,
      z = ~ umap_3,
      marker = list(
        size = size,
        color = ~ tmp.expr,
        colorscale = color.pal,
        showscale = TRUE),
      text = ~ paste("Cluster: ", Cluster)) %>%
      add_markers() %>%
      layout(xaxis = ax, yaxis = ax)
    
    htmlwidgets::saveWidget(as_widget(plot), file)
  } else {
    pData(cds)$tmp.expr = log2(cds@assayData$exprs[gene_short_name,] / pData(cds)$Size_Factor + 1)
    
    ax <- list(
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
    )


    plot = plot_ly(
      pData(cds)[ids,],
      x = ~ umap_1,
      y = ~ umap_2,
      z = ~ umap_3,
      marker = list(
        size = size,
        color = ~ tmp.expr,
        colorscale = color.pal,
        showscale = TRUE),
      text = ~ paste("Cluster: ", Cluster, "<br>log2 norm expr: ", tmp.expr)) %>%
      add_markers() %>%
      layout(xaxis = ax, yaxis = ax)
    
    pData(cds)$tmp.expr = NULL
    
    htmlwidgets::saveWidget(as_widget(plot), file)
  }
}

# get a table of mean expression, fraction expressing, and specifity of all genes per cell type
aggregated_expr_data = function (cds, group_cells_by = "cell_type_all"){

cell_group_df <- data.frame(row.names = row.names(colData(cds)), 
                            cell_id = row.names(colData(cds)))
cell_group_df$cell_group <- colData(cds)[, group_cells_by]
cell_group_df$cell_group <- as.character(cell_group_df$cell_group)

cluster_binary_exprs = as.matrix(aggregate_gene_expression(cds, 
                  cell_group_df = cell_group_df, norm_method = "binary"))

cluster_fraction_expressing_table = tibble::rownames_to_column(as.data.frame(cluster_binary_exprs))
cluster_fraction_expressing_table = tidyr::gather(cluster_fraction_expressing_table, 
                                                  "cell_group", "fraction_expressing", -rowname)

cluster_mean_exprs = as.matrix(aggregate_gene_expression(cds, 
                              cell_group_df = cell_group_df, norm_method = "size_only"))

cluster_expr_table = tibble::rownames_to_column(as.data.frame(cluster_mean_exprs))
cluster_expr_table = tidyr::gather(cluster_expr_table, "cell_group", 
                                   "mean_expression", -rowname)

cluster_fraction_expressing_table$mean_expression = cluster_expr_table$mean_expression

cluster_spec_mat = monocle3:::specificity_matrix(cluster_mean_exprs, cores = 4)
cluster_spec_table = tibble::rownames_to_column(as.data.frame(cluster_spec_mat))
cluster_spec_table = tidyr::gather(cluster_spec_table, "cell_group", 
                                   "specificity", -rowname)

cluster_fraction_expressing_table$specificity = cluster_spec_table$specificity
cluster_fraction_expressing_table = cluster_fraction_expressing_table %>% 
      dplyr::rename("gene_id" = rowname) %>% 
      dplyr::left_join(rowData(cds) %>% 
                      as.data.frame() %>%
                      dplyr::select("gene_id" = "id", gene_short_name), 
              by = "gene_id") %>% 
      dplyr::select(cell_group, gene_id, gene_short_name, everything())

return(cluster_fraction_expressing_table)
}


### perform DEA between two sets of cells (from JP)
two.set.differential.gene.test.m3 = function(cds, set.1.filter, set.2.filter, formal = F, cores = 1, thresh = 1) {
  message(paste("# of cells in set 1:", sum(set.1.filter)))
  message(paste("# of cells in set 2:", sum(set.2.filter)))
  
  s1.cds = cds[, set.1.filter]
  s2.cds = cds[, set.2.filter]
  
  s1.norm.expr = sweep(assay(s1.cds), 2, colData(s1.cds)$Size_Factor, "/")
  s2.norm.expr = sweep(assay(s2.cds), 2, colData(s2.cds)$Size_Factor, "/")
  
  s1.tpm = apply(s1.norm.expr, 1, sum)
  s1.tpm = s1.tpm / sum(s1.tpm) * 1000000
  s2.tpm = apply(s2.norm.expr, 1, sum)
  s2.tpm = s2.tpm / sum(s2.tpm) * 1000000
  
  s1.n.umi = apply(assay(s1.cds), 1, sum)
  s2.n.umi = apply(assay(s2.cds), 1, sum)
  
  higher.expr = ifelse(s1.tpm > s2.tpm, "Set 1", "Set 2")
  
  s1.ratio = s1.tpm / (s2.tpm + 1)
  s2.ratio = s2.tpm / (s1.tpm + 1)
  log2.ratio = ifelse(
    s1.tpm == 0 & s2.tpm == 0, 0, ifelse(
      higher.expr == "Set 1", log2(s1.ratio), log2(s2.ratio)))
  
  s1.n.expr = apply(assay(s1.cds), 1, function(x) sum(x >= thresh))
  s2.n.expr = apply(assay(s2.cds), 1, function(x) sum(x >= thresh))
  
  s1.precision = s1.n.expr / (s1.n.expr + s2.n.expr)
  s1.recall = s1.n.expr / ncol(s1.cds)
  s2.precision = s2.n.expr / (s1.n.expr + s2.n.expr)
  s2.recall = s2.n.expr / ncol(s2.cds)
  
  precision = ifelse(higher.expr == "Set 1", s1.precision, s2.precision)
  recall = ifelse(higher.expr == "Set 1", s1.recall, s2.recall)
  
  f.score = 2 * precision * recall / (precision + recall)
  
  res = data.frame(
    gene = rowData(cds)$id,
    symbol = rowData(cds)$gene_short_name,
    set.1.n.umi = s1.n.umi,
    set.2.n.umi = s2.n.umi,
    set.1.tpm = s1.tpm,
    set.2.tpm = s2.tpm,
    higher.expr = higher.expr,
    log2.ratio = log2.ratio,
    precision = precision,
    recall = recall,
    f.score = f.score
  ) %>% arrange(-f.score)
  
  
  
  if (formal) {
    colData(cds)$tmp = ifelse(set.1.filter, 1, ifelse(set.2.filter, 2, NA))
    
    cds.subset = cds[, set.1.filter | set.2.filter]
    # cds.subset = estimateSizeFactors(cds.subset)
    # cds.subset = estimateDispersions(cds.subset)
    cds.subset = detect_genes(cds.subset, 0.1)
    
    expressed.genes = subset(rowData(cds.subset), num_cells_expressed >= 5)[, 1]
    message(paste(length(expressed.genes), "genes expressed in at least 5 cells across both sets"))
    message("Computing differential expression p-values")
    
    DEG = monocle::differentialGeneTest(cds.subset[expressed.genes,],
                               fullModelFormulaStr = "~ tmp", cores = cores)
      
    print(head(DEG))
    
    res = inner_join(res,
                     DEG %>% select(gene = id,
                                    p.val = pval, q.val = qval),
                     by = "gene") %>% arrange(q.val)
    
    colData(cds)$tmp = NULL
  }
  
  return(res)
}


get.norm.expr.matrix = function(cds) {
    mat = assay(cds)
    mat@x = mat@x / rep.int(colData(cds)$Size_Factor, diff(mat@p))
    return(mat)
}

append_umap_coordinates = function(cds, umap_3D = F){
  
  if (!umap_3D){
    colData(cds)$umap1 = reducedDim(x = cds,
                                    type = "UMAP")[,1]
    colData(cds)$umap2 = reducedDim(x = cds,
                                    type = "UMAP")[,2]
  }
  if (umap_3D){
    colData(cds)$umap3d_1 = reducedDim(x = cds,
                                       type = "UMAP")[,1]
    colData(cds)$umap3d_2 = reducedDim(x = cds,
                                       type = "UMAP")[,2]
    colData(cds)$umap3d_3 = reducedDim(x = cds,
                                       type = "UMAP")[,3]
  }
  return(cds)
}
