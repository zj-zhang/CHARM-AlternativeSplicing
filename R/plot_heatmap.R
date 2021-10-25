# A wrapper for use in jupyter notebook for plotting heatmaps
source("R/heatmap.3.R")

fill_row_with_mean = function(x)
{
    filled_na = 0
    for(i in 1:nrow(x)) {
        idx = which(is.na(x[i,]))
        filled_na = filled_na + length(idx)
        row_mean = mean(x[i,], na.rm=T)
        row_mean = ifelse(is.na(row_mean), 0, row_mean)
        x[i, idx] = row_mean
    }
    cat("filled_na=",filled_na,"\n")
    return(x)
}


plot_heatmap = function(plot_data, clab, main="heatmap.3", rescale=TRUE) {
    plot_data = as.matrix(plot_data)
    if(rescale) {
	y = t(scale(t(plot_data)))
    } else {
	y = plot_data
    }
    print(dim(y))
    y = fill_row_with_mean(y)

    match_index = match(colnames(y), rownames(clab))
    clab = as.matrix(clab[match_index,])
    has_var_row = which(apply(y, 1, var)>0)
    print(length(has_var_row))
    y = y[has_var_row, ]

    hr <- hclust(as.dist(1-cor(t(y), method="pearson", use="complete")), method="average")
    hc <- hclust(as.dist(1-cor(y, method="pearson", use="complete")), method="average")
    #hr <- hclust(dist(y, method="euclidean"), method="average")
    #hc <- hclust(dist(t(y), method="euclidean"), method="average")

    palette <- colorRampPalette(c("yellow3","white","darkblue"))

    breaks = c(seq(-2,-1,length=100), seq(-0.99, 0.99,length=100), seq(1, 2,length=100))

    print(dim(y))

    heatmap.3(y, Rowv = as.dendrogram(hr), Colv = as.dendrogram(hc), dendrogram = "both", col = palette,
      ColSideColors = clab, 
      key = TRUE,
      breaks=breaks,
      main=main
    )
}

