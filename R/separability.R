##################### Separability ####################

#' separability_pairwise
#'
#' Calculate separability for every sample pair. The higher the more separable.
#'
#' @param tsne_coord
#' A dataframe of the two dimension t-SNE coordinates of each cell in the combined sample.
#' @param group
#' A vector of group assignment for each cell, with the same order as the tsne_coord.
#' @param sample_label
#' A vector of sample labels for each cell, with the same order as the tsne_coord.
#' @param k
#' K-nearest neighbours used to calculate distance. DEFAULT is 5.
#' @return A matrix of separability for each sample pair (column) and each group (row).
#' @export


separability_pairwise <- function(tsne_coord, group, sample_label, k = 5)
{## separability_pairwise will call function separability_by_group.
	message("calculate separability for each sample pair")
	sample_pair <- combn(levels(sample_label), 2)
	colnames(sample_pair) <- paste0(sample_pair[1, ], '.vs.', sample_pair[2, ])
	res <- apply(sample_pair, 2, function(x){
		print(x)
		ind <- sample_label %in% x
		tsne_coord_sub <- tsne_coord[ind, ]
		group_sub <- group[ind]
		sample_label_sub <- sample_label[ind]
		sepa <- separability_by_group(tsne_coord_sub, group_sub, sample_label_sub, k = k)
		return(sepa)
	})
	return(res)
}

#' separability_by_group
#'
#' Calculate separability for each group in a pair of samples. Internal function called by separability_pairwise.
#'
#' @param tsne_coord
#' A dataframe of the two dimension t-SNE coordinates of each cell in the combined sample.
#' @param group
#' A vector of group assignment for each cell, with the same order as the tsne_coord.
#' @param sample_label
#' A vector of sample labels for each cell, with the same order as the tsne_coord.
#' @param k
#' K-nearest neighbours used to calculate distance. DEFAULT is 5.
#' @return A vector of separability for each group.
#' @export

separability_by_group <- function(tsne_coord, group, sample_label, k = 5)
{ ## separability_by_group will call function separability.
    p <- sapply(levels(group), function(i){
        m <- tsne_coord[group == i, ]
		if (nrow(m) == 0) avg_diff <- NA else
		{
			class_label <- sample_label[group == i]
			avg_diff <- separability(m, class_label, k = k)
		}
    })
    p <- p/(max(tsne_coord[, 1]) - min(tsne_coord[, 1]))*100
    round(p, 2)
}

#' separability
#'
#' Calculate separability for labeled data. Internal function called by separability_by_group.
#'
#' @param coord
#' A dataframe of the two dimension t-SNE coordinates of each cell in the combined sample.
#' @param class_label
#' A vector of sample labels for each cell, with the same order as the tsne_coord.
#' @param k
#' K-nearest neighbours used to calculate distance. DEFAULT is 5.
#' @return A single value of separability.
#' @export

separability <- function(coord, class_label, k = 5)
{ ## separability will call function inter_dist and inna_dist.
    if(length(unique(class_label)) > 1)
    {
        d <- dist(coord)
        dls <- split(as.data.frame(as.matrix(d)), class_label)
        sample1_dls <- split(as.data.frame(t(dls[[1]])), class_label)
        d1 <- sample1_dls[[1]] ## distance matrix between cells in sample1
        d21 <- sample1_dls[[2]] ## distance matrix from cells in sample2 to cells in sample1
        sample2_dls <- split(as.data.frame(t(dls[[2]])), class_label)
        d2 <- sample2_dls[[2]] ## distance matrix between cells in sample2
        d12 <- sample2_dls[[1]] ## distance matrix from cells in sample1 to cells in sample2
        diff1 <- inter_dist(d12, k)-inna_dist(d1, k)
        diff2 <- inter_dist(d21, k)-inna_dist(d2, k)
        avg_diff <- mean(c(diff1, diff2)) ## average two samples
    } else avg_diff <- Inf
    return(avg_diff)
}

#' inter_dist
#'
#' Calculate the inter-sample distance. Internal function called by separability.
#'
#' @param x
#' A distance matrix of cells in one sample to the cells in another sample.
#' @param k
#' K-nearest neighbours used to calculate distance. DEFAULT is 5.
#' @return A single value of distance.

inter_dist <- function(x, k)
{
    ## knn_mean <- function(x, k) {mean(sort(x)[1:k])}
    ## mean(apply(x, 1, knn_mean, k))
    knn_median <- function(x, k){ median(sort(x)[1:k]) } ## take median of K distance of a cell
    median(apply(x, 1, knn_median, k)) ## take median of distance of all cells
}

#' inna_dist
#'
#' Calculate the inna-sample distance. Internal function called by separability.
#'
#' @param x
#' A distance matrix of cells in one sample to the cells in the same sample.
#' @param k
#' K-nearest neighbours used to calculate distance.
#' @return A single value of distance.

inna_dist <- function(x, k)
{
    ## knn_mean <- function(x, k) {mean(sort(x)[2:k+1])}
    ## mean(apply(x, 1, knn_mean, k))
    knn_median <- function(x, k){ median(sort(x)[2:k+1]) } ## take median of K distance of a cell
    median(apply(x, 1, knn_median, k)) ## take median of distance of all cells
}

