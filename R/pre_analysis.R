##################### Pre-analysis ####################

#' make_single_obj
#'
#' A warper of Seurat function to generate Seurat object and marker genes for single sample.
#' @import Seurat
#' @param data_dir
#' Directory with 10X genomics single cell data or full path of expression table.
#' @param is.10X
#' The input data is 10X genomics format or not. DEFAULT is TRUE.
#' @param output
#' The output directory to save the marker genes.
#' @return A Seurat object.
#' @export

make_single_obj <- function(data_dir, is.10X = TRUE, output)
{

	if (is.10X) da <- Read10X(data_dir) else da <- read.table(data_dir, sep = "\t")
	obj <- CreateSeuratObject(raw.data = da, min.cells = 3, min.genes = 200)
	obj <- NormalizeData(obj)
	obj <- FindVariableGenes(obj)
	obj <- ScaleData(obj, vars.to.regress = c("nUMI"))
	obj <- RunPCA(obj, do.print = FALSE)
	obj <- RunTSNE(obj, dims.use = 1:10)
	obj <- FindClusters(obj, dims.use = 1:10, resolution = 0.6)
	markers <- FindAllMarkers(obj, only.pos = TRUE)
	write.csv(markers, file = paste0(output, '.markers.csv'))
	return(obj)
}

#' make_comb_obj
#'
#' A warper of Seurat function to generate Seurat object for combined sample from single sample.
#' @import Seurat
#' @param data_dirList
#' A list of directory with 10X genomics single cell data or full path of expression table.
#' @param is.10X
#' The input data is 10X genomics format or not. DEFAULT is TRUE.
#' @param comb_delim
#' The delimiter used in the cell names in the combined object to connect sample name and cell name in individual sample. DEFAULT is '-'.
#' @return A Seurat object.
#' @export

make_comb_obj <- function(data_dirList, is.10X = TRUE, comb_delim = '-')
{

	## rename cell names in each dataset and then combine
	n <- names(data_dirList)
	if (is.null(n)) {
		names(data_dirList) <- paste0('s', 1:length(data_dirList))
		warning("The names(data_dirList) is NULL. Sample names are assigned as '", paste(names(data_dirList), collapse = ' '), "'")
		}
	da_list <- lapply(n, function(x){
		if (is.10X) da <- Read10X(data_dirList[[x]]) else da <- read.table(data_dirList[[x]], sep = "\t")
		colnames(da) <- paste0(x, comb_delim, colnames(da))
		return(da)
		})
	comb_da <- do.call(cbind, da_list)

	obj <- CreateSeuratObject(raw.data = comb_da, min.cells = 3, min.genes = 200)
	obj <- NormalizeData(obj)
	obj <- FindVariableGenes(obj)
	obj <- ScaleData(obj, vars.to.regress = c("nUMI"))
	obj <- RunPCA(obj, do.print = FALSE)
	obj <- RunTSNE(obj, dims.use = 1:10)
	return(obj)
}
