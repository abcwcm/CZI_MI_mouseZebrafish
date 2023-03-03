library(magrittr);library(data.table)

#' Add reduced dimensions from one SCE to another
#' 
#' @details This is particularly useful when subsets of a data set have been
#' re-clustered and re-redDim'ed. The function here takes care of merging
#' the data sets and adding NAs for missing cells.
#' 
#' @param sce1 SCE object to which a new redDim should be added
#' @param sce2 SCE object from which the new redDim will be taken
#' @param rd2 which new redDim of sce2 should be added?
#' @return SCE object with added redDim
#' @examples
#' scall <- add_UMAP(scall, sc9, rd2 = "UMAP", prefix = "cluster9_10only")
add_UMAP <- function(sce1 = scall, sce2 = sc9, rd1 = "UMAP", rd2 = "UMAP",
    prefix = "sc9"){
    rdold <- reducedDim(sce1, rd1) %>% as.data.table(., keep.rownames=TRUE) 
    rdnew <- reducedDim(sce2, rd2) %>% as.data.table(., keep.rownames=TRUE)
    rdnew <- rdnew[rdold, on = "rn"]
    
    newx <- paste(prefix, "x", sep = ".")
    newy <- paste(prefix, "y", sep = ".")
    
    setnames(rdnew, names(rdnew), c("Row.names", newx, newy,"old.x","old.y"))
    rdnew <- as.data.frame(rdnew)
    rownames(rdnew) <- rdnew$Row.names
    rdnew <- rdnew[, c(newx,newy)]
    
    scout <- sce1
    reducedDim(scout, paste(rd2, prefix, sep = "_")) <- rdnew[colnames(scout),]
    return(scout)
}
