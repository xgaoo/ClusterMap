#######################################################
####################  ClusterMap  #####################
#######################################################

#################### Master function ##################

#' cluster_map
#'
#' A master function to perform the full workflow of ClusterMap.
#'
#' @import ggplot2
#' @import pheatmap
#' @import ape
#' @import Seurat
#' @import circlize
#'
#' @importFrom grDevices col2rgb dev.off hcl dev.off pdf png rgb
#' @importFrom graphics plot
#' @importFrom stats dist hclust median setNames
#' @importFrom utils combn read.csv read.table write.csv
#'
#' @param marker_file_list
#' A list of csv files with names. Each file is a marker gene table for a sample. The columns named as 'cluster' and 'gene' are required.
#' @param edge_cutoff
#' The edge length cutoff to decide the sub-nodes to merge or not. DEFAULT is 0.1.
#' @param output
#' The output directory to save the matching results.
#' @param cell_num_list
#' A list of vector of cell numbers for each group and each sample.
#' @param single_obj_list
#' A list of Seurat object for each sample, with the same list names as the list names of marker_file_list.
#' @param comb_obj
#' A Seurat object for the combined sample. Cells in different samples are labelled by the sample names with the comb_delim. The sample names should be the same as the list names of marker_file_list.
#' @param comb_delim
#' The delimiter used in the cell names in the combined object to connect sample name and cell name in individual sample. DEFAULT is '-'.
#' @param k
#' K-nearest neighbours used to calculate distance. DEFAULT is 5.
#' @param reduction
#' Select the reduction of "tsne", "umap", or "pca" that used for the recolor image.
#' @return A dataframe of the matching results. Heatmap of marker genes, the corresponding dendrogram, circos plot and recolored t-SNE plots will be saved into files.
#' @export


cluster_map <- function(marker_file_list, edge_cutoff = 0.1, output, cell_num_list = NULL, single_obj_list = NULL, comb_obj = NULL, comb_delim = '-', k = 5, seurat_version = 3, reduction="tsne")
{
	circos.clear()
	## Version check for comb delim
	if(seurat_version == 3){
		comb_delim = '_'
	}
	## match sub groups
	mapRes <- cluster_map_by_marker(marker_file_list, cutoff = edge_cutoff, output = output)

	## pull out cell_num_list if single Seurat object list is provided.
	if (!is.null(single_obj_list))
	{
		if (all(names(marker_file_list) == names(single_obj_list)) == FALSE | is.null(names(marker_file_list)) | is.null(names(single_obj_list)))
			stop("names(marker_file_list) doesn't match names(single_obj_list).")

        if (single_obj_list[[1]]@version > 3){
		cell_num_list <- lapply(single_obj_list,
                            function(obj){
                                    summary(Idents(obj))
                                         })
        }

        else if(single_obj_list[[1]]@version < 3){
		cell_num_list <- lapply(single_obj_list, function(obj){
                                    summary(obj@ident)})
	    }
    }
	## make circos plot and add cell percentage if cell_num_list is provided or single Seurat object list is provided.
	if (!is.null(cell_num_list))
	{
		if (all(names(marker_file_list) == names(cell_num_list)) == FALSE | is.null(names(marker_file_list)) | is.null(names(cell_num_list)))
			stop("names(marker_file_list) doesn't match names(cell_num_list).")
		circos_map(mapRes, cell_num_list, output)
		mapRes <- add_perc(mapRes, cell_num_list)
	}

	## Recolor reduction plot for each sample if single Seurat object list is provided.
	if (!is.null(single_obj_list))
	{
		sample_names <- names(single_obj_list)
		new_group_list <- lapply(sample_names, function(n){
			da <- structure(as.vector(mapRes[, n]), names = mapRes$regroup)
			recolor_s(da, single_obj_list[[n]], n, reduction=reduction)
		})
		return(new_group_list)
		names(new_group_list) <- names(single_obj_list)

		## Recolor reduction plot for combined sample and calculate separability if combined Seurat object is provided.
		if (!is.null(comb_obj))
		{
			sample_label <- as.factor(sub(paste0(comb_delim, '.*'), '', rownames(comb_obj@meta.data)))
			if (all(sort(levels(sample_label)) == sort(names(new_group_list))) == FALSE)
				{
				 message("Sample label in comb_obj: ")
				 print(levels(sample_label))
				 message("names(single_obj_list): ")
				 print(names(single_obj_list))
				 stop("Sample label in comb_obj doesn't match names(new_group_list) or names(single_obj_list).")
				 }
			
			new_group_list$comb <- recolor_comb(comb_obj, new_group_list, output, comb_delim, reduction=reduction)

			coord <- as.data.frame(comb_obj@reductions[[reduction]]@cell.embeddings)
			sepa <- separability_pairwise(coord, group = new_group_list$comb, sample_label, k = k)
			colnames(sepa) <- paste0(colnames(sepa), '_separability')

			mapRes <- cbind(mapRes, sepa)
		}
		saveRDS(new_group_list, file = paste0(output, '.new.group.list.RDS'))
	}
	write.csv(mapRes, file = paste0(output, '.results.csv'))
	return(mapRes)
}
