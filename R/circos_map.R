##################### Circos plot #####################

#' circos_map
#'
#' Plot Circos plot for the matching results.
#'
#' @import circlize
#' @param mapRes
#' A dataframe of the output of function cluster_map_by_marker.
#' @param cell_num_list
#' A list of vector of cell numbers for each group and each sample.
#' @param output
#' The output directory to save the plot.
#' @param color_cord
#' A vector of colors for the cord of circos plot. DEFAULT is NULL. Pre defined internal color will be used.
#' @param color_sample
#' A vector of colors for the sample sectors in the circos plot. DEFAULT is NULL. Pre defined internal color will be used.
#' @param width
#' The width of circos plot by inch.
#' @param height
#' The height of circos plot by inch.
#' @return circos plot will be save.
#' @export


circos_map <- function(mapRes, cell_num_list, output, color_cord = NULL, color_sample = NULL, width=7, height=7)
{ ## circos_map will call function plot_circos, gg_color_hue and makeTransparent.
	message("circos plot")

	sample_name <- names(cell_num_list)
	if (all(sample_name %in% colnames(mapRes)) == FALSE)
		stop("names(marker_file_list) or samples in mapRes doesn't match names(cell_num_list).")
	## pairwise link
	combs <- combn(1:length(sample_name), 2)
	mapRes_samples <- mapRes[, sample_name]
	colnames(mapRes_samples) <- NA
	pair <- do.call(rbind, apply(combs, 2, function(x) mapRes_samples[, x]))
	pair <- cbind(pair, mapRes[, c('regroup', 'similarity')])
	pair <- pair[pair[, 1] != '' & pair[, 2] != '' & !is.na(pair[, 1]) & !is.na(pair[, 2]), ]
	colnames(pair)[1:2] <- c('v1', 'v2')
	## split merged clusters with ;
	pair_list <- lapply(1:nrow(pair), function(i)
	{
		y <- pair[i, ]
		v1 <- unlist(strsplit(as.vector(y$v1), ';'))
		v2 <- unlist(strsplit(as.vector(y$v2), ';'))
		v <- c()
		for (x in v1) v <- rbind(v, cbind(x, v2))
		v <- as.data.frame(v, stringsAsFactors = FALSE)
		v$similarity <- y$similarity
		v$regroup <- y$regroup
		return(v)
	})
	pair <- do.call(rbind, pair_list)
	colnames(pair)[1:2] <- c('v1', 'v2')
	pair$v1 <- as.vector(pair$v1)
	pair$v2 <- as.vector(pair$v2)
	## cell number percentage
	cell_perc_list <- lapply(cell_num_list, function(x) round(x/sum(x), 2))
	if (is.null(color_sample)) col_sample <- rep(c("#f9865c", "#84e281", "#74d2f7", "#e083fc", "#ffbf66", "#6682ff"), 10) else col_sample <- color_sample
	if (is.null(color_cord)) col_cord <- gg_color_hue(nrow(mapRes)) else col_cord <- color_cord
	## circos plot
	png(paste0(output, '.circos.png'), width=width/7*480, height=height/7*480)
		plot_circos(cell_perc_list, pair, mapRes, col_cord, col_sample)
	dev.off()
	pdf(paste0(output, '.circos.pdf'), width=width, height=height)
		plot_circos(cell_perc_list, pair, mapRes, col_cord, col_sample)
	dev.off()
}


plot_circos <- function(cell_perc_list, pair, mapRes, col_cord, col_sample)
{
	temp=cell_perc_list
	names(temp)=paste0(names(temp),"__")
	cell_perc <- unlist(temp)
	names(cell_perc) <- sub('__.', '__', names(cell_perc))
	cell_perc[cell_perc < 0.01] <- 0.01 ## too small to plot
	fa <- factor(names(cell_perc), levels = unique(names(cell_perc)))
	## initialize
	gaps <- lapply(cell_perc_list, function(x) c(rep(1, length(x)-1), 8))
	circos.par(gap.after = unlist(gaps), start.degree = -3, cell.padding = c(0, 0, 0, 0))
	circos.initialize(fa, xlim = cbind(rep(0, length(cell_perc)), cell_perc))
	## plot sample sectors
	circos.track(ylim = c(0, 1), track.height = uh(5, "mm"), bg.border = NA)
	for (n in names(cell_perc_list)) highlight.sector(paste0(n, "__", names(cell_perc_list[[n]])), track.index = 1,
		col = col_sample[match(n, names(cell_perc_list))], padding = c(0, 0, 0.3, 0), text = n, cex = 1.5, text.col = "black", niceFacing = TRUE)
	## plot sub-group sectors
	circos.track(fa, ylim = c(0, 1), panel.fun = function(x, y)
		{
		circos.text(CELL_META$xcenter, CELL_META$ylim[1], sub('.*__', '', CELL_META$sector.index),
		adj = c(0.3, -2), niceFacing = TRUE)
		}, bg.col = 'black', bg.border = NA, track.height = 0.05, track.margin = c(0, 0.1)
	)
	## plot cord
	for(i in 1:nrow(pair))
	{
		x <- pair[i, ]
		col <- makeTransparent(col_cord[x$regroup], round(x$similarity*100))
		circos.link(x$v1, c(0, cell_perc[x$v1]), x$v2, c(0, cell_perc[x$v2]), col = col, border = NA, h.ratio = 0.5)
	}
	circos.clear()
}


gg_color_hue <- function(n)
{
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}


makeTransparent <- function(someColor, alpha = 100)
{
    col <- col2rgb(someColor)
    rgb(red = col[1], green = col[2], blue = col[3], alpha = alpha, maxColorValue = 255)
  }

########## Add cell percentage to result table ########

#' circos_map
#'
#' Plot Circos plot for the matching results.
#'
#' @import circlize
#' @param mapRes
#' A dataframe of the output of function cluster_map_by_marker.
#' @param cell_num_list
#' A list of vector of cell numbers for each group and each sample.
#' @return A dataframe of the matching results with cell percentage column added.
#' @export

add_perc <- function(mapRes, cell_num_list)
{
	sample_name <- names(cell_num_list)
	if (all(sample_name %in% colnames(mapRes)) == FALSE)
		stop("names(marker_file_list) or samples in mapRes doesn't match names(cell_num_list).")
	cell_perc_list <- lapply(cell_num_list, function(x) round(x/sum(x), 2))
	temp=cell_perc_list
	names(temp)=paste0(names(temp),"__")
	cell_perc <- unlist(temp)
	names(cell_perc) <- sub('__.', '__', names(cell_perc))
	## add to mapRes
	res_sub <- mapRes[, sample_name]
	res_perc <- apply(res_sub, 1:2, function(x)
	{
		paste(cell_perc[unlist(strsplit(x, ';'))], collapse = ';')
	})
	colnames(res_perc) <- paste0(colnames(res_perc), '_cell_perc')
	res <- cbind(mapRes, res_perc)
	return(res)
}
