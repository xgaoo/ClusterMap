################# Clustering and tree cut #############

#' cluster_map_by_marker
#'
#' Match groups by marker genes and decompose into new groups by purity cut.
#'
#' @import pheatmap
#' @import ape
#' @param marker_file_list
#' A list of csv files. Each file is a marker gene table for a sample. The columns named as 'cluster' and 'gene' are required.
#' @param cutoff
#' The edge length cutoff to decide the sub-nodes to merge or not. DEFAULT is 0.1.
#' @param output
#' The output directory to save the matching results.
#' @return A dataframe of the matching results. Heatmap of marker genes and the dendrogram will be saved into files.
#' @export

cluster_map_by_marker <- function(marker_file_list, cutoff = 0.1, output)
{ ## cluster_map_by_marker will call function purity_cut.
	message("match sub-clusters")
	## get marker table
	if (is.null(names(marker_file_list)))
	{
		names(marker_file_list) <- paste0('s', 1:length(marker_file_list))
		warning("The names(marker_file_list) is empty. Sample names are assigned as '", paste(names(marker_file_list), collapse = ' '), "'" )
	}
	markerList <- lapply(marker_file_list, read.csv, as.is = T)
	markerList <- lapply(names(markerList), function(n)
	{
		x <- markerList[[n]]
		x$cluster <- paste0(n, '__', x$cluster)
		return(x)
	})
	markers <- do.call(rbind, markerList)
	## clustering
	da <- table(markers[,c("cluster","gene")])
	d <- dist(da, method = 'binary')
	hc <- hclust(d, method = 'average')
	png(paste0(output, '.hcluster.png'), width = 480*length(marker_file_list)) ## save dendrogram png
		plot(hc)
	dev.off()
	pdf(paste0(output, '.hcluster.pdf'), width = 7*length(marker_file_list)) ## save dendrogram pdf
		plot(hc)
	dev.off()
	## save heatmap png
	png(paste0(output, '.heatmap.hcluster.png'), height = 480*length(marker_file_list))
	ph <- pheatmap(da, scale = 'none', clustering_method = 'average', color = c('#e3f8f9', '#fc2807'),
		show_rownames = T, show_colnames = F, clustering_distance_rows = 'binary',
		legend_breaks = c(0, 1), legend_labels = c(0, 1))
	dev.off()
	ph <- pheatmap(da, scale = 'none', clustering_method = 'average', color = c('#e3f8f9', '#fc2807'),
		show_rownames = T, show_colnames = F, clustering_distance_rows = 'binary',
		legend_breaks = c(0, 1), legend_labels = c(0, 1),
		filename = paste0(output, '.heatmap.hcluster.pdf'), height = 7*length(marker_file_list))
	## tree cut
	res <- purity_cut(hc, cutoff)
	##write.csv(res, file = paste0(output, '.cluster.map.csv'))
	return(res)
}

#' purity_cut
#'
#' Cut hierarchical clustering dendrogram by edge length and purity of the nodes.
#'
#' @import ape
#' @param hcluster
#' A hclust object.
#' @param cutoff
#' The edge length cutoff to decide the sub-nodes to merge or not. DEFAULT is 0.1.
#' @return A dataframe of the matching results.
#' @export

purity_cut <- function(hcluster, cutoff = 0.1)
{
	hcp <- as.phylo(hcluster)
	png('ape.tree.png')
		plot(hcp, edge.width = 2, label.offset = 0.1)
		nodelabels()
		tiplabels()
	dev.off()
	## get edge info
	tree <- cbind(hcp$edge, round(hcp$edge.length, 3))
	colnames(tree) <- c('high_node', 'low_node', 'edge_length')
	tree <- as.data.frame(tree, as.is = T)
	tree$too_long <- (tree$edge_length > cutoff/2) ## edge_length is height/2
	## get all the offspring of each node
	offs_nodeList <- lapply(tree$low_node, function(x) {
    if (x <= length(hcp$tip.label)) return(hcp$tip.label[x]) else
		return(extract.clade(hcp, x)$tip.label)})
	names(offs_nodeList) <- tree$low_node
	## get sample list that the nodes belong to
	sampleList <- lapply(offs_nodeList, function(x) unique(sub('__.*$', '', x)))
	tree$low_node_offs_sample <-  unlist(lapply(sampleList, length))
	## check if keep the node
	nodeList <- split(tree, tree$high_node)
	n_single <- unlist(lapply(nodeList, function(x) sum(x$low_node_offs_sample == 1)))
	n_too_long <- unlist(lapply(nodeList, function(x) sum(x$too_long)))
	no_include <- unlist(lapply(nodeList, function(x)
	{
		low_node <- as.character(x$low_node)
		sampleList_sub <- sampleList[low_node]
		include <- (all(sampleList_sub[[1]] %in% sampleList_sub[[2]]) | all(sampleList_sub[[2]] %in% sampleList_sub[[1]]))
		return(!include)
	}))
	keep <- (n_single == 2 | (n_single == 1 & n_too_long < 2)| (n_single == 0 & n_too_long == 0 & no_include))
	names(keep) <- names(nodeList)
	## check if any offspring node is cut
	node_merge_order <- rev(unique(tree$high_node))
	for (nd in node_merge_order)
	{
		nd <- as.character(nd)
		offs_nodes <- nodeList[[nd]]$low_node
		keep_sub <- keep[as.character(offs_nodes)]
		keep_sub[is.na(keep_sub)] <- TRUE
		tmp <- (sum(keep_sub) == 2 & keep[nd])
		keep[nd] <- tmp
	}
	if (all(keep == TRUE))
        {
                message("No matched groups.")
                return(NULL)
        }

	## get rid of lower duplicated nodes
	res <- offs_nodeList[names(keep)[keep == T]]
	for (x in rev(names(res)))
	{
		if (any(res[[x]] %in% unlist(res[names(res) != x]))) res = res[names(res) != x]
	}
	## get singletons
	singles <- setdiff(hcp$tip.label, unlist(res))
	res <- c(res, as.list(singles))
	## reform output
	lev <- sort(unique(sub('__.*$', '', hcp$tip.label)))
	res_reform <- do.call(rbind, lapply(res, function(x)
	{
		tmp <- split(x, f = factor(sub('__.*$', '', x), levels = lev))
		unlist(lapply(tmp, paste, collapse = ';'))
	}))
	rownames(res_reform) <- 1:nrow(res_reform)
	res_reform[res_reform == ''] <- NA
	res_reform <- as.data.frame(res_reform, stringsAsFactors = FALSE)
	## add similarity
	res_reform$similarity <- round(1 - branching.times(hcp)[names(res)]*2, 2)
	## add group
	res_reform <- res_reform[order(-res_reform$similarity), ]
	res_reform$regroup <- 1:nrow(res_reform)
	rownames(res_reform) <-  res_reform$regroup
	return(res_reform)
}
