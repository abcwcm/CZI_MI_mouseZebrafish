#' Extract a subset of cells from a composite SCE file
#' 
#' @details The only point of this function is to harmonize the 
#' resulting SCE files to enable the generation of a list of SCEs
#' that are ready for MNN-integration, i.e. they have counts() etc.
#' @param insce SCE object
#' @return SCE object
#' @examples
#' scel <- lapply(unique(sce$Sample), function(x) individualize_samples(sce[, sce$Sample == x]))
#' names(scel) <- unique(sce$Sample)
individualize_samples <- function(insce, pf = "old",
    shared_colData = c("cell","Sample","Barcode","Tissue","condition","post.surgery")){
    
    if(! ( all( shared_colData %in% names(colData(insce)) ) )){
        stop("Check the colData columns you've indicated. They don't all seem to be present.")
    }
    
    outsce <- insce
    colData(outsce) <- colData(outsce)[, shared_colData]
    logcounts(outsce) <- NULL
    if("UMAP" %in% reducedDimNames(outsce)){
        reducedDim(outsce, paste0("UMAP_", pf))<- reducedDim(outsce, "UMAP")
        reducedDim(outsce, "UMAP") <- NULL
    }
    if("PCA_corr" %in% reducedDimNames(outsce)){
        reducedDim(outsce, paste0("PCA_corr_", pf)) <- reducedDim(outsce, "PCA_corr")
        reducedDim(outsce, "PCA_corr") <- NULL
    }
    # rownames(outsce) <- scater::uniquifyFeatureNames(rowData(outsce)$ID, rowData(outsce)$Symbol)
    rownames(outsce) <- rowData(outsce)$ID
    return(outsce)
}


#' Perform MNN-integration, clustering and UMAP'ing
#' @param list_of_sce
#' @param hvg_n number of highly variable genes. Default: 2000.
#' @param shared_colData names of the columns to be extracted that are present in
#' all SCEs of \code{list_of_sce}.
#' @param cluster_ks integer values for k = neighborhood setting for the clustering
#' @param fn filename prefix
#' @param save_hvg whether to save the objects related to the highly variable genes.
#' @param ignore_genes vector of gene names to be ignored for the MNN'ing
#' Default: TRUE
MNN_correct <- function(list_of_sce, hvg_n = 2000,
    shared_colData = c("cell","Sample","Barcode","Tissue","condition","post.surgery"),
    fn="allCells_ZF_sampleIntegration", 
    cluster_ks = c(100,200),
    save_hvg = TRUE, ignore_genes=NULL){
        
        if(!(all(unlist(lapply(list_of_sce, function(x) all(shared_colData %in% names(colData(x)))))))){
            stop("Check the colData columns you've indicated. They don't seem to be present in all of the SCEs of list_of_sce")
        }
        if(!all(unlist(lapply(list_of_sce, function(x) all(c("ID","Symbol") %in% names(rowData(x))))))){
            stop(print(" 'ID' and 'Symbol' should be part of the rowData"))
        }
        
        universe <- Reduce(intersect, lapply(list_of_sce, rownames))
        list_of_sce2 <- lapply(list_of_sce, "[", i=universe)
        # generate logcounts
        normed.sce <- do.call(multiBatchNorm, list_of_sce2) # returns a list
        # Identifying a set of HVGs using stats from all batches, using logcounts
        all.dec <- lapply(normed.sce, modelGeneVar)
        combined.dec <- do.call(combineVar, all.dec)
        combined.hvg <- getTopHVGs(combined.dec, n=hvg_n) 
        
        if(save_hvg == TRUE){
            save(all.dec, combined.hvg, combined.hvg, file = paste0(fn, "_HVGs.rda"))
        }
        
        if(!is.null(ignore_genes)){
            if(any(ignore_genes %in% combined.hvg)){
                rm_genes <- ignore_genes[ignore_genes %in% combined.hvg]
                combined.hvg <- combined.hvg[!combined.hvg %in% rm_genes]
                }else if(!any(ignore_genes %in% universe)){
                    print("None of the genes to be ignored were part of original data. Check the gene name format.")
                }
        }
        
        # Merge with MNN ----------------------------
        ## prep
        combined <- noCorrect(normed.sce)
        assayNames(combined) <- "logcounts"
        combined$Sample <- combined$batch
        set.seed(1010100)
        ## progressively merge cells from each sample in each batch until all cells 
        ## are mapped onto a common coordinate space
        print("Performing MNN")
        multiout <- fastMNN(combined, batch=combined$Sample, subset.row=combined.hvg)
        # Renaming metadata fields for easier communication later.
        multiout$Sample <- multiout$batch
        ## UMAP-----------------------------------
        print("Performing UMAP")
        set.seed(10101010)
        multiout <- runUMAP(multiout, dimred="corrected")
        saveRDS(multiout,file = paste0(fn, "_MNNmerged.rds"))
        
        ## CLUSTERING -----------------------------
        print(paste0("Performing clustering for ", paste("k =", cluster_ks, collapse = ", " )))
        for(i in cluster_ks){
            g <- buildSNNGraph(multiout, use.dimred="corrected", k = i)
            clusters <- igraph::cluster_louvain(g)
            multiout[[paste0("cluster_k",i)]] <- factor(clusters$membership)
        }

        ## combine
        universe <- Reduce(intersect, lapply(list_of_sce, rownames))
        list_of_sce <- lapply(list_of_sce, "[", i=universe)
        comb.mat <- lapply(list_of_sce, counts) %>% do.call(cbind, .)
        
        ## colData 
        
        cd <- lapply(list_of_sce, function(x) colData(x)[, shared_colData]) %>% 
            do.call(rbind, .)
        
        ## add clusters from multiout to combined SCE 
        cd2 <- as.data.frame(colData(multiout))
        cd2$cell <- rownames(cd2)
        cd2 <- cd2[, c("cell",grep("cluster_k", names(cd2), value=TRUE))]
        newcd <- merge(cd, cd2, by = "cell") 
        rownames(newcd) <- newcd$cell
        newcd <- newcd[colnames(comb.mat),]
        
        ### rowData
        rd <- rowData(list_of_sce[[1]])[, c("ID","Symbol")]
        rd <- rd[rownames(comb.mat),]
        
        out.sce <- SingleCellExperiment(
            assays = list(counts = comb.mat),
            colData = newcd, rowData = rd)
        
        ## add redDims from the merged data set
        rdu <- reducedDim(multiout, "UMAP") 
        reducedDim(out.sce, "UMAP") <- rdu[colnames(out.sce),]
        reducedDim(out.sce, "PCA_corr") <- reducedDim(multiout, "corrected")
        
        ## add log-counts
        qckclst <- quickCluster(out.sce, method = "igraph", min.mean = 0.1)
        out.sce <- computeSumFactors(out.sce, min.mean=0.1, cluster = qckclst)
        out.sce <- scater::logNormCounts(out.sce)
        
        ## minor tweaks
        is.mito <- grepl("mt-", ignore.case = TRUE, rowData(out.sce)$Symbol)
        out.sce <- scater::addPerCellQC(out.sce, subsets=list(mito=is.mito))
        out.sce$log10_total_features <- log10(out.sce$detected)
        
        return(out.sce)
}
