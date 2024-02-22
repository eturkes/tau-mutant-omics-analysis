#    This file is part of tau-mutant-omics-analysis.
#    Copyright (C) 2023-2024  Emir Turkes, Naoto Watamura, UK DRI at UCL
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Emir Turkes can be contacted at emir.turkes@eturkes.com

# This file holds common functions and methods.

#' ggplot2 function providing custom aesthetics and automatic placement of categorical labels.
#' For continuous data, a colorbar is implemented.
#'
#' @param data SingleCellExperiment or Seurat object.
#' @param x,y Dimensionality reduction coordinates.
#' @param color Column metadata to color points by.
#' @param type \code{"cat"} is categorical, \code{"cont"} is continuous, \code{"NULL"} is generic.
#' @examples
#' red_dim_plot(data = sce, x = "tsne1", y = "tsne2", color = "cluster", type = "cat")
#' red_dim_plot(data = seurat, x = "umap1", y = "umap2", color = "nUMI", type = "cont")
#'
red_dim_plot <- function(data, x, y, color, type = NULL) {

  if ((class(data))[1] == "SingleCellExperiment") {
    gg_df <- data.frame(colData(data)[ , c(x, y, color)])
  } else if ((class(data))[1] == "Seurat") {
    gg_df <- data.frame(data[[x]], data[[y]], data[[color]])
  }
  rownames(gg_df) <- NULL
  gg_df[[color]] <- factor(gg_df[[color]])

  gg <- ggplot(gg_df, aes_string(x, y, col = color)) +
    geom_point(alpha = 0.35, stroke = 0.05, shape = 21, aes_string(fill = color)) +
    theme_classic() +
    theme(
      legend.position = "right", plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1)))

  if (is.null(type)) {
    return(gg)

  } else if (type == "cat") {
    label_df <- gg_df %>% group_by_at(color) %>% summarise_at(vars(x:y), median)
    label_df <- cbind(label_df[[1]], label_df)
    names(label_df) <- c("label", color, x, y)
    gg <- gg + geom_label_repel(
      data = label_df, max.overlaps = Inf, aes(label = label), show.legend = FALSE
    )

  } else if (type == "cont") {
    if ((class(data))[1] == "SingleCellExperiment") {
      gg_df <- data.frame(colData(data)[ , c(x, y, color)])
    } else if ((class(data))[1] == "Seurat") {
      gg_df <- data.frame(data[[x]], data[[y]], data[[color]])
    }
    rownames(gg_df) <- NULL

    gg <- ggplot(gg_df, aes_string(x, y)) +
      geom_point(alpha = 0.35, stroke = 0.05, aes_string(color = color)) +
      theme_classic() +
      theme(
        legend.position = "right", plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()
      ) +
      scale_color_viridis()
  }
  gg
}

#' Adds download buttons and horizontal scrolling to \code{"DT::datatable"}.
#'
#' @param dt A data.table object.
#' @examples
#' datatable_download(dt = data_table)
#'
datatable_download <- function(dt) {

  datatable(
    dt,
    list(
      scrollX = TRUE, dom = "Blfrtip",
      buttons = list(
        "copy", "print",
        list(extend = "collection", buttons = c("csv", "excel", "pdf"), text = "Download")
      )
    ),
    extensions = "Buttons"
  )
}

#' Adds download buttons, horizontal scrolling, exponential values to \code{"DT::datatable"}.
#'
#' @param dt A data.table object.
#' @examples
#' datatable_download_exp(dt = data_table)
#'
datatable_download_exp <- function(dt) {

  datatable(
    dt,
    list(
      scrollX = TRUE,
      dom = "Blfrtip",
      buttons = list(
        "copy", "print",
        list(extend = "collection", buttons = c("csv", "excel", "pdf"), text = "Download")
      ),
      rowCallback = JS(
        "function(row, data) {",
        "for (i = 1; i < data.length; i++) {",
        "if (data[i]>=1000 | data[i]<1000) {",
        "$('td:eq('+i+')', row).html(data[i].toExponential(2));}}}"
      )
    ),
    extensions = "Buttons"
  )
}

#' Convert human to mouse gene names.
#' Adapted from:
#' https://rjbioinformatics.com/2016/10/14/converting-mouse-to-human-gene-names-with-biomart-package/
#'
#' @param genes A vector of human genes.
#' @examples
#' human_to_mouse_genes(genes = gene_list)
#'
human_to_mouse_genes <- function(genes) {

  human <- useMart("ensembl", "hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", "mmusculus_gene_ensembl")

  new_genes <- getLDS(
    attributes = "ensembl_gene_id", filters = "ensembl_gene_id", values = genes,
    mart = human, attributesL = "external_gene_name", martL = mouse
  )
  new_genes <- unique(new_genes[ , 2])
}

#' Pipeline for normalization, dimensionality reduction, and clustering of post-QC scRNA-seq data.
#'
#' @param seurat Post-QC Seurat object.
#' @param cache_dir Directory to save post-processed Seurat object.
#' @param sub_name Subset level for naming of cache object.
#' @param protocol Vector with the following elements in this order: \code{"human"} or
#' \code{"mouse"}. \code{"droplet"} or \code{"smart-seq"}. \code{"single-cell"} or
#' \code{"single-nuc"}. \code{"umis"} or \code{"reads"}.
#' @param vars_to_regress Vector of nuisance variables for sctransform to regress out.
#' @param parallel_override See function \code{"parallel_plan"}.
#' @param cc Logical, whether to perform cell-cycle scoring.
#' @param res_divider, What to divide number of cells by to arrive at clustering resolution.
#' @param conserve_memory, Logical, whether to use conserve.memory in sctransform.
#' @param min_cells, Numerical, minimum number of cells to retain a gene in sctransform.
#' @examples
#' cluster_pipeline(
#'   seurat = seurat, cache_dir = cache_dir, sub_name = "neuronal",
#'   protocol = protocol, vars_to_regress = "mito_percent", parallel_override = NULL,
#'   cc = FALSE, res_divider = 1000, conserve.memory = TRUE
#' )
cluster_pipeline <- function(
    seurat, cache_dir, sub_name, protocol, vars_to_regress,
    parallel_override, cc = TRUE, res_divider = 3000, conserve_memory = FALSE, min_cells = NULL
) {

  rds <- file.path(cache_dir, paste0("seurat_", sub_name, ".rds"))
  if (file.exists(rds)) {
    seurat <- readRDS(rds)
    return(seurat)
  } else {

    if (protocol[4] == "umis") {
      # Run sctransform.
      # Note that this function produces many iterations of the following benign warning:
      # Warning in theta.ml(y = y, mu = fit$fitted): iteration limit reached
      # ---------------------------------------------------------------------------------
      parallel_plan(seurat, parallel_override)
      if (conserve_memory == TRUE & !is.null(min_cells)) {
      seurat <- suppressWarnings(
        SCTransform(
            seurat, vars.to.regress = vars_to_regress, vst.flavor = "v2",
            conserve.memory = TRUE, min_cells = min_cells, verbose = FALSE
          )
        )
      } else if (conserve_memory == TRUE) {
        seurat <- suppressWarnings(
          SCTransform(
            seurat, vars.to.regress = vars_to_regress, vst.flavor = "v2",
            conserve.memory = TRUE, verbose = FALSE
          )
        )
      } else if (conserve_memory == FALSE & !is.null(min_cells)) {
        seurat <- suppressWarnings(
          SCTransform(
            seurat, vars.to.regress = vars_to_regress, vst.flavor = "v2", min_cells = min_cells, verbose = FALSE
          )
        )
      } else if (conserve_memory == FALSE) {
        seurat <- suppressWarnings(
          SCTransform(
            seurat, vars.to.regress = vars_to_regress, vst.flavor = "v2", verbose = FALSE
          )
        )
      }
      # ---------------------------------------------------------------------------------

      # Perform PCA.
      # ------------
      seurat <- RunPCA(seurat, verbose = FALSE)
      add_df <- data.frame(Embeddings(seurat)[ , 1:2])
      names(add_df) <- paste0("pca", seq(ncol(add_df)))
      seurat$pca1 <- add_df$pca1
      seurat$pca2 <- add_df$pca2
      reduction <- "pca"
      dims <- 1:30 # Dimensions for downstream computations.
      # ------------

    } else if (protocol[4] == "reads") {
      sce <- as.SingleCellExperiment(seurat) # ZINB-WaVE can directly take an SCE but not Seurat.
      logcounts(sce) <- NULL

      # Use top 1,000 variable features for downstream computations.
      # ------------------------------------------------------------
      seurat <- NormalizeData(seurat, verbose = FALSE)
      seurat <- FindVariableFeatures(seurat, nfeatures = 1000, verbose = FALSE)
      sce <- sce[rownames(sce) %in% VariableFeatures(seurat), ]
      # ------------------------------------------------------------

      # Get ENSEMBL annotations for more accurate retrieval of gene length and GC content.
      # ----------------------------------------------------------------------------------
      if (protocol[1] == "human") {
        dataset <- "hsapiens_gene_ensembl"
      } else if (protocol[1] == "mouse") {
        dataset <- "mmusculus_gene_ensembl"
      }
      mart <- useEnsembl("ensembl", dataset)
      attributes <- c("external_gene_name", "ensembl_gene_id")
      gene_anno <- getBM(attributes, "external_gene_name", rownames(sce), mart)
      # ----------------------------------------------------------------------------------

      # For gene symbols with multiple ENSEMBL IDs, duplicate the gene symbol to have an identical
      # row for each ENSEMBL ID.
      # ------------------------------------------------------------------------------------------
      dup <- gene_anno[duplicated(gene_anno$external_gene_name), ]
      for (i in 1:dim(dup)[1]) {
        for (j in 1:dim(gene_anno)[1]) {
          if (dup$ensembl_gene_id[i] == gene_anno$ensembl_gene_id[j]) {
            gene_anno$external_gene_name[j] <- paste0(gene_anno$external_gene_name[j], "-dup")
          }
        }
      }
      sce <- sce[rownames(sce) %in% gene_anno$external_gene_name, ]
      new_mat <- counts(sce)
      for (i in 1:dim(dup)[1]) {
        for (j in 1:dim(sce)[1]) {
          if (dup$external_gene_name[i] == rownames(sce)[j]) {
            new_row <- t(counts(sce)[j, ])
            rownames(new_row) <- paste0(rownames(counts(sce))[j], "-dup")
            new_mat <- rbind(new_mat, new_row)
          }
        }
      }
      gene_anno <- gene_anno[gene_anno$external_gene_name %in% rownames(new_mat), ]
      gene_anno <- gene_anno[order(match(gene_anno$external_gene_name, rownames(new_mat))), ]
      rownames(new_mat) <- gene_anno$ensembl_gene_id
      sce <- SingleCellExperiment(list(counts = new_mat), colData = colData(sce))
      # ------------------------------------------------------------------------------------------

      # Get gene lengths and GC content.
      # --------------------------------
      row_data <- data.frame(getGeneLengthAndGCContent(rownames(sce), dataset, "biomart"))
      rowData(sce) <- data.frame(gc_content = row_data$gc, length = row_data$length)
      # --------------------------------

      # Run ZINB-WaVE.
      # TODO: Fix hardcoding of `vars_to_regress`.
      # ------------------------------------------
      counts(sce) <- as.matrix(counts(sce)) # ZINB-WaVE is incompatible with dgCMatrix.
      sce <- zinbwave(
        sce, paste0("~ ", vars_to_regress[1], " + ", vars_to_regress[2]), "~ gc_content + length",
        10, epsilon = 1e12, BPPARAM = MulticoreParam()
      )

      zinb <- data.frame(reducedDim(sce, "zinbwave"))
      seurat$zinb1 <- zinb$W1
      seurat$zinb2 <- zinb$W2
      zinb <- as.matrix(zinb)
      seurat[["zinb"]] <- CreateDimReducObject(zinb, key = "W", assay = DefaultAssay(seurat))
      reduction <- "zinb"
      dims <- 1:10 # Dimensions for downstream computations.
      # ------------------------------------------
    }

    # Perform cell cycle scoring.
    # ---------------------------
    if (cc == TRUE) {
      if (protocol[1] == "human") {
        s_genes <- cc.genes.updated.2019$s.genes
        g2m_genes <- cc.genes.updated.2019$g2m.genes
      } else if (protocol[1] == "mouse") {
        s_genes <- human_to_mouse_genes(cc.genes.updated.2019$s.genes)
        g2m_genes <- human_to_mouse_genes(cc.genes.updated.2019$g2m.genes)
      }
      seurat <- CellCycleScoring(seurat, s_genes, g2m_genes)
      seurat$cc_diff <- seurat$S.Score - seurat$G2M.Score # Combined proliferating cell signal.
    }
    # ---------------------------

    # Perform UMAP reduction.
    # -----------------------
    seurat <- RunUMAP(seurat, dims, reduction, verbose = FALSE)
    add_df <- data.frame(Embeddings(seurat, "umap"))
    names(add_df) <- paste0("umap", seq(ncol(add_df)))
    seurat$umap1 <- add_df$umap1
    seurat$umap2 <- add_df$umap2
    # -----------------------

    # Perform Louvain clustering.
    # ---------------------------
    resolution <- ncol(seurat) / res_divider
    seurat <- FindNeighbors(seurat, reduction, dims, verbose = FALSE)
    seurat <- FindClusters(seurat, resolution = resolution, verbose = FALSE)
    # ---------------------------

    saveRDS(seurat, rds)
  }
  seurat
}

#' Set the \code{"plan"} for \code{"future"} based on free memory and object size with the option to
#' override.
#'
#' @param object Object to check if \code{"future.globals.maxSize"} large enough to parallelize.
#' @param parallel_override \code{"NULL"} to calculate plan decision, \code{0} for sequential, a
#' non-zero integer for multiprocess and to set \code{"future.globals.maxSize"}.
#' @examples
#' parallel_plan(object = seurat, parallel_override = 5368709120)
#'
parallel_plan <- function(object, parallel_override = NULL) {

  if (is.null(parallel_override)) {
    # Get free memory.
    # ----------------
    gc()
    mem <- as.numeric(unlist(strsplit(system("free -b", TRUE)[2], " "))[15])
    # ----------------

    # Distribute free memory (minus 10 GiB) across available cores.
    # -------------------------------------------------------------
    mem <- mem - 10 * 1024 ^ 3
    mem <- mem / availableCores()
    # -------------------------------------------------------------

    # Enable parallelization only if `object` can fit in `future.globals.maxSize` (plus 1 Gib).
    # -----------------------------------------------------------------------------------------
    if (mem > object.size(object) + 1 * 1024 ^ 3) {
      plan("multisession")
      options(future.globals.maxSize = mem)
    } else {
      plan("sequential")
    }
    # -----------------------------------------------------------------------------------------

  } else if (parallel_override == 0) {
    plan("sequential")

  } else {
    plan("multisession")
    options(future.globals.maxSize = parallel_override)
  }
}

computeGeneSetsOverlapMax <- function(gSets, uniqGenes, min.sz=1, max.sz=Inf) {
  ## gSetsMembershipMatrix should be a (genes x gene-sets) incidence matrix

  gSetsMembershipMatrix <- incidence(gSets)
  gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, colnames(gSetsMembershipMatrix) %in% uniqGenes])

  lenGsets <- colSums(gSetsMembershipMatrix)

  szFilterMask <- lenGsets >= max(1, min.sz) & lenGsets <= max.sz
  if (!any(szFilterMask))
    stop("No gene set meets the minimum and maximum size filter\n")

  gSetsMembershipMatrix <- gSetsMembershipMatrix[, szFilterMask]
  lenGsets <- lenGsets[szFilterMask]

  totalGsets <- ncol(gSetsMembershipMatrix)

  M <- t(gSetsMembershipMatrix) %*% gSetsMembershipMatrix

  M1 <- matrix(lenGsets, nrow=totalGsets, ncol=totalGsets,
               dimnames=list(colnames(gSetsMembershipMatrix), colnames(gSetsMembershipMatrix)))
  M2 <- t(M1)
  M.max <- matrix(0, nrow=totalGsets, ncol=totalGsets)
  M.max[M1 > M2] <- M1[M1 > M2]
  M.max[M2 >= M1] <- M2[M2 >= M1]
  overlapMatrix <- M / M.max

  return (overlapMatrix)
}

word_cloud = function(x, width = NULL, dict = c()){
  t.rW = c(c("the", "and", "that", "process",
           "binding", "cell", "which", "activity",
           "from", "protein", "any", "complex",
           "with", "response", "proteins", "are",
           "organism", "enzyme", "movement", "change",
           "results", "membrane", "state", "receptor",
           "stimulus", "production", "result", "leading",
           "into", "etc", "form", "cells",
           "its", "expression", "specific", "terms",
           "within", "for", "involved", "rna",
           "none", "group", "catalysis", "series",
           "transport", "two", "found", "ion",
           "reaction", "formation", "signaling", "reactions",
           "assembly", "acid", "structure", "such",
           "chemical", "nuclear", "small", "between",
           "pathways", "resulting", "cellular", "pathway",
           "molecule", "also", "molecular", "contains",
           "other", "part", "amino", "out",
           "binds", "including", "some", "outer",
           "directed", "signal", "intracellular", "signals",
           "complexes", "phosphate", "associated", "during",
           "membranes", "type", "via", "factor", "mature",
           "ring", "large", "metabolic",
           "plasma", "target", "activation", "composed",
           "both", "more", "acids", "arrangement",
           "involving", "this", "body", "gene",
           "one", "hydrolysis", "lumen", "domain",
           "transfer", "surface", "acts", "action",
           "ribonucleoprotein", "has", "level", "can",
           "through", "region", "cytoplasm",
           "molecules", "system", "residues", "growth",
           "localization", "their", "ligand", "ions",
           "together", "where", "whose", "generated",
           "proteincoupled", "processing", "most", "specialized",
           "pore", "progression", "time", "single",
           "formed", "constituent", "contain", "transporter",
           "active", "nucleus", "components", "several",
           "have", "over", "another",
           "location", "side", "consisting", "hydrogen",
           "often", "central", "forms", "downstream",
           "modulates", "maintenance", "role", "electron",
           "may", "site", "catalyzes", "subunit",
           "containing", "after", "organization", "when",
           "disassembly", "reduces", "maintained", "receptors",
           "acceptor", "regulation", "family", "cycle",
           "structures", "substances", "include", "inner",
           "initial", "snrnp", "compound", "ends",
           "initiated", "while", "outcome", "those",
           "removal", "class", "smooth", "function",
           "attached", "regulator", "compounds", "acting",
           "means", "across", "features", "tissue",
           "andor", "free", "enzymes", "activated",
           "regulatory", "increase", "sequence", "they",
           "either", "includes", "parts", "enables",
           "not", "many", "tissues", "position",
           "related", "acquires", "usually", "donor",
           "these", "subunits", "messenger", "initiation",
           "volume", "occurs", "carried", "mediated",
           "component", "well", "but", "bonding",
           "increases", "set", "transported", "release",
           "agent", "concentration", "known", "ending",
           "precursor", "modification", "play", "core",
           "export", "substance", "biosynthetic", "different",
           "enclosed", "end", "synthesized", "located",
           "groups", "species", "light",
           "portion", "example", "synthesis", "subsequent",
           "assembled", "released", "certain", "unique",
           "defined", "metabolism", "head", "heptameric",
           "producted", "envelope", "elongation", "organized",
           "cytosol", "each", "products", "factors",
           "esterified", "interaction", "distal", "four",
           "addition", "remain", "produced", "involves",
           "donors", "consists", "similar", "chain",
           "second", "cytosolic", "relatively", "exposure",
           "apparatus", "unspecialized", "information", "internal",
           "distinct", "chains", "dense", "substrate",
           "physiological", "functional", "nervous", "consequence",
           "water", "interacting", "bound", "organisms",
           "direct", "particle", "diphosphate", "first",
           "primary", "three", "hydrophobic", "secretion",
           "secretory", "sites", "residue", "mechanism",
           "will", "cone", "population", "external",
           "alpha", "following", "member", "due",
           "levels", "covalent", "humans", "number",
           "contributes", "biosynthesis", "carbon", "cytoplasmic",
           "required", "attachment", "fusion",
           "beta", "expansion", "processes", "all",
           "derived", "occuring", "bonds", "present",
           "appearance", "derivative", "mammals", "soluble",
           "rise", "bond", "called", "only",
           "hydroxyl", "transduction", "organs", "using",
           "composition", "occuring", "high", "gamma",
           "template", "step", "eukaryotic", "naromatic",
           "double", "substrates", "surrounding", "functions",
           "anion", "having", "host", "forming",
           "length", "hydroxy", "highly", "domains",
           "combining", "polymerase", "initiate", "typically",
           "establishment", "adhesion", "interactions", "organelle",
           "shortening", "presence", "uatac", "characterized",
           "prevents", "aromatic", "monoterpenoid", "dismantling",
           "conversion", "derivable", "diketones", "terpenoids",
           "necessary", "pyruvate", "quinone", "rearrangement",
           "mass", "csmp", "size", "errorfree",
           "postreplication", "prr", "psyp", "heat",
           "isomerase", "phosphoric", "symbiont", "thyroxine",
           "doublet", "line", "total", "effect",
           "cap", "strand", "targets", "metal",
           "amount", "retained", "drugs", "killing",
           "transports", "bip", "proteincontaining", "recruitment",
           "superoxide", "tfiia", "selenium", "commonly",
           "amide", "largely", "induce", "triiodothyronine",
           "polymers", "copper", "affect", "stimulation",
           "own", "unwinding", "branch",
           "duplex", "linkages", "added", "continues",
           "proportion", "separate", "longitudinal", "controls",
           "elastic", "proximity", "pass", "respectively",
           "rupture", "alcohol", "reformation", "meats",
           "tyrosinebased", "express", "evidence", "microorganisms",
           "codon", "aptype", "basic",
           "associates", "electrondense", "spherical", "close",
           "turn", "lie", "rather", "programmed",
           "rather", "break", "endogenous", "lens",
           "regrowth", "main", "however", "singlestrand",
           "mouth", "amides", "produce", "carboxylic",
           "interacts", "main", "watersoluble", "pole",
           "mainly", "upstream", "fructose", "stops",
           "main", "produce", "rather", "modifications",
           "case", "characteristics", "predominantly", "ligation",
           "asymmetry", "extends", "especially", "doublestrand",
           "branches", "glycosidic", "aqueous", "branching",
           "takes", "primarily", "tubular", "membranous",
           "staining", "ribbon", "bending", "perceived",
           "maternal", "interior", "upregulation", "should",
           "interconnected", "centers", "relative", "principal",
           "multiprotein", "existing", "movements", "prior",
           "entering", "promotes", "consist", "spontaneous",
           "conjugating", "circulating", "corresponding", "glutamylcysteinylglycine",
           "polymerization", "hydrolase", "extending", "allows",
           "next", "common", "transforming", "combination",
           "nmultimeric", "dinucleotide", "cues", "linked",
           "life", "adenine", "guidance", "positively",
           "mus", "microscopy", "replacement", "diffuse",
           "cactivating", "multimeric", "responses", "prenyl",
           "integrate", "musculus", "contiguous", "guanine",
           "shaping", "recognizes", "isoprenoids", "klinked",
           "duct", "white", "fourth", "chronic",
           "lead", "guanylnucleotide", "map", "shaft",
           "ventral", "received", "posterior", "preinitiation",
           "received", "pic", "contribute", "hybrid",
           "clustering", "aldehyde", "uptake", "diffusion",
           "long", "particular", "condition", "noncovalent",
           "condensation", "like", "messengers", "still",
           "very", "homologous", "lower", "leftright",
           "lining", "mechanical", "intramolecular", "lineages",
           "glycerol", "attain", "amine", "stimulate",
           "described", "left", "daughter", "sorting",
           "aliphatic", "way", "direction", "density",
           "stacks", "arrival", "colonystimulating", "attractive",
           "gravity", "independent", "alditols", "successive",
           "networks", "dependent", "permitting", "ionic",
           "repulsive", "broken", "flat", "order",
           "coding", "heterotrimeric", "poles", "block",
           "sigma", "gyrus", "secreted", "capacity",
           "iron", "contributing", "tolerance", "acidic",
           "homology", "limiting", "adult", "secreted",
           "liquid", "astral", "inositol", "atpdependent",
           "diverse", "phosphatidyl", "sensor", "approximately",
           "leucine", "foods", "normal", "solutes",
           "transition", "betaalanine", "ligands", "normal",
           "range", "originates", "collateral", "product",
           "emerge", "perinuclear", "internalization", "belonging",
           "transmitting", "purpose", "residing", "transmitting",
           "condensed", "lineage", "sometimes", "disulfide",
           "transmit", "effects", "induced", "measured",
           "releasing", "intake", "down", "less",
           "ammonia", "differ", "nearly", "loop",
           "threonine", "just", "interact", "driven",
           "bundle", "compact", "fluids", "associate",
           "allowing", "carries", "complexed", "context",
           "udpnacetylgalactosamine", "naturally", "destined",
           "cleaves", "proper", "brings", "bipolar",
           "diacylglycerol", "extension", "protrusion", "comprises",
           "critical", "plant", "plane", "above",
           "causes", "poisonous", "grow", "overall",
           "destruction", "there", "until", "force",
           "take", "mediate", "assembles", "migrating",
           "rounded", "enable", "chemomechanical", "playing",
           "acyclic", "chco", "provides", "given",
           "retains", "coupled", "freely", "spliced",
           "geometry", "canal", "chnh", "dimer",
           "alcoholic", "matter", "links", "sulfhydryl",
           "radial", "delivered", "recruit", "subcomplex",
           "away", "dimers", "spindle", "promoter",
           "transcribed", "compartment", "integrin", "tripeptide",
           "new", "surrounds", "mannose", "classii",
           "specialization", "recycling", "opposite", "presentation",
           "segment", "toward", "kinases",
           "resting", "phospholipid", "maintaining", "alteration",
           "does", "organizing", "ability", "glycerophospholipids",
           "accumulation", "electrons", "tubes", "motile",
           "sheath", "making", "expresses", "abundance",
           "self", "whole", "cisternae",
           "partner", "stability", "periodicity", "represent",
           "intrinsic", "execution", "exist", "restricted",
           "moving", "ribonucleoside", "monounsaturated", "aminoacyltrna",
           "bilayers", "cadmium", "polymer",
           "committed", "betalinked", "gdp", "providing",
           "factorbeta", "marker", "element", "separated",
           "polysaccharides", "variable", "underlying", "asymmetric",
           "point", "tethering", "seen", "fbox",
           "accessory", "agents", "alternatives", "basement",
           "distinctive", "lamellae", "actions", "inorganic",
           "giant", "dynamic", "reach", "regulated",
           "synthetic", "increasing", "animals", "doublestranded",
           "being", "requires", "tract", "initially",
           "lack", "name", "determination", "nodes",
           "rich", "anisotropic", "perform", "upper",
           "ribose", "five", "third", "bridging",
           "temperature", "then", "bridge", "pair",
           "homo", "considered", "channels", "near",
           "store", "carrying", "filled", "recognize",
           "sequences", "perception", "attaches", "rest",
           "positions", "lserine", "assemble", "clusters",
           "sorted", "ceramide", "continuous", "elongated",
           "adherens", "potentials", "atom", "amounts",
           "reactive", "initiator", "protons", "ternary",
           "inducing", "deacetylase", "alternative", "dioxygenase",
           "gradient", "atoms", "intact", "full",
           "constituents", "periphery", "sets", "soma",
           "originating", "roles", "aminoacid", "deposited",
           "mineralized", "connection", "combined", "substituted",
           "sequential", "neurological", "perpetuation", "crawling",
           "original", "solution", "faces", "produces",
           "carboxyterminal", "tight", "carbonyl", "primitive",
           "although", "cdc", "kda", "dnatemplated",
           "transferase", "precursors", "mixture", "opening",
           "creation", "additional", "somatic", "energyindependent",
           "effectors", "polysaccharide", "protontransporting", "lobe",
           "tube", "plate", "dorsal", "axes",
           "patterns", "anteriorposterior", "particles", "comprise",
           "recruits", "essential", "bridges", "elements",
           "thick", "regulators", "hkme", "peroxidase",
           "pentose", "border", "enter", "irreversibly",
           "observed", "bisphosphate", "arranged", "cis",
           "electrochemical", "oligosaccharyltransferase", "types", "once",
           "innervate", "clear", "proceeds", "rough",
           "nucleolar", "conveyed", "spans", "dissociation",
           "come", "granular", "depolymerization", "completion",
           "exonucleolytic", "anchored", "ammonium", "creating",
           "minus", "antiparallel", "sides", "share",
           "them", "dioxide", "establish", "correct",
           "changed", "sulfurcontaining", "antigenindependent", "leaflet",
           "proteindna", "removes", "functioning", "thioester",
           "thought", "cortical", "processed", "final",
           "nonidentical", "facilitated", "guanylyl", "terminal",
           "cavity", "internally", "whereby", "nterminal",
           "sheet", "alphamannosyltransferase", "anhydride", "covers",
           "exit", "lacks", "udpnacetyldglucosamine", "tubule",
           "carbons", "selenocysteine", "air", "fashion",
           "novo", "stage", "extended", "sensation",
           "crucial", "individuals", "deprivation", "areas",
           "carboxylate", "units", "zone", "membranebounded",
           "ctype", "acidity", "aminohydroxypropanoic", "nitrate",
           "halves", "shu", "arms", "acquire",
           "modulate", "agonist", "gammadelta", "elsewhere",
           "among", "properties", "bind", "denantiomer",
           "alkaline", "serine", "specification", "association",
           "stable", "labyrinth", "prevented", "dserine",
           "tie", "basicity", "oligosaccharide", "orientation",
           "give", "linoleic", "phenotypes", "measure",
           "antagonist", "capable", "acceptors", "udpsugar",
           "eventually", "mononucleate", "right", "wide",
           "been", "alone", "respect", "higher",
           "under", "place", "solute", "linkage",
           "spoke", "acetoacetate", "rhythm", "selfrenewing",
           "acylphosphatas", "plays", "major", "switching",
           "widely", "hydroxybutyrate", "biochemically", "acylphosphatase",
           "cholate", "sphinganine", "synthase", "instead",
           "galactosylceramide", "transiently", "hydroxybutanoate", "ceramides",
           "choloylcoa", "reduce", "prenylated", "abdominal",
           "hours", "integrity", "galactosylceramidesulfate", "rhydroxybutanoate",
           "autonomic", "regularity", "palate", "deoxygalactose",
           "dorsalventral", "frequent", "affected", "detection",
           "intermediates", "every", "nacetyllactosamine", "dead",
           "biologically", "elasticity", "recur", "identity",
           "activities", "front", "classical", "repellence",
           "runs", "prenylatedprotein", "opens", "resistance",
           "phosphatidylglycerol", "phosphooligonucleotides", "comprised", "polya",
           "spr", "twostage", "galbetaglcnacbetan", "polynacetyllactosamine",
           "aminocarboxypropyltransferase", "proceeding", "basal", "specifically",
           "monomers", "polarized", "udpgalactose", "globular",
           "negative", "nitric", "blue", "flavan",
           "shock", "protrusions", "joining", "overhang",
           "remaining", "nine", "deiminase", "removing",
           "according", "loss", "oxidase", "phenol",
           "become", "fluid", "fuses", "superfamily",
           "adpribosyltransferase", "nain", "proteinlysine", "naout",
           "early", "sodiumsulfate", "sulfatein", "sulfateout",
           "symporter", "barrier", "referred", "crystalline",
           "potential", "aum", "thin", "induces",
           "acylchain", "multiple", "elimination", "scallopshaped",
           "nacetylaspartylglutamate", "antimicrobial", "adaptor",
           "origin", "bulk", "superficial", "thickened",
           "nacetyllaspartate", "key", "umbrella", "nacetyllaspartatelglutamate",
           "collects", "gets", "beside", "uroplakin",
           "excretion", "get", "around", "proteinaceous",
           "selective", "alphadgalactosylbetadgalactosylbetanacetyldglucosaminylr", "sulfide", "alphagalactosyltransferase",
           "smaller", "shup", "triggered", "acceleration",
           "afferent", "angular", "oxide", "electromagnetic",
           "methionine", "wavelength", "flavanoids", "hphosphothreonine",
           "nuclei", "flavones", "hthreonine", "flavonols",
           "remove", "beds", "intercellular", "activating",
           "khdckhdcl", "modulation", "noncovalently", "commitment",
           "selectively", "monoester", "optimum", "nlrp",
           "subcortical", "tle", "incorporated", "alphabeta",
           "blocks", "arabinonucleosides", "constitutively", "held",
           "dsn", "arabinopentose", "expressed", "helps",
           "adaptors", "arabinose", "mis", "biology",
           "mismind", "cascades", "darabinose", "nsl",
           "pmf", "larabinose", "unit", "about",
           "glycocalyx", "outermost", "rnas", "aminogroups",
           "adhesive", "apicalbasolateral", "helix", "deactivation",
           "cushion", "generate", "formaldehyde", "cranially",
           "retrograde", "bleaching", "methanethiol", "ganglionic",
           "hydrolytic", "prevent", "yellow", "augment",
           "integrated", "projects", "edge", "lascorbic",
           "stalk", "whereas", "positive", "surrounded",
           "lines", "rapidly", "cerebral", "anterior",
           "reductase", "bodies", "characteristic", "neural",
           "striated", "arginine", "droplet", "space",
           "animal", "satellite", "acidification", "stem",
           "origins", "nucleoside", "trigger", "center",
           "rectification", "carriermediated", "monomeric", "diphosphatase",
           "receive", "ndc", "constitute", "laid",
           "byproducts", "wraps", "heavy", "lateral",
           "made", "inward", "trimer", "subapical",
           "provide", "doublets", "oxygen", "occurring",
           "passage", "adjacent", "contact", "modified",
           "throughout", "reduced", "important", "regeneration",
           "serinethreonine", "excretory", "followed", "catalyzed",
           "filtration", "converted", "bcdk", "aminoacylation",
           "coupling", "derivatives", "synthetase", "transferred",
           "increased", "piece", "connected", "tip",
           "compartments", "arylamine", "black", "numbers",
           "brown", "aminocarboxypropyluridine", "nacetylarylamine", "methylthioadenosine",
           "nacetyltransferase", "trnauridine", "act", "human",
           "adherent", "multicellular", "classes", "ultimately",
           "variety", "migrate", "situated", "oaminoacyltrna",
           "neuronal", "sulfate", "transesterification", "classi",
           "carboxyl", "genes", "stimuli", "phosphoadenosine",
           "phosphosulfate", "acyl", "extraocular", "pattern",
           "stimulates", "cascade", "conjugated", "nerve",
           "fucose", "guanosine", "heterodimeric", "phenylpropanoid",
           "morphologically", "phosphorus", "nonantigenic", "network"),
           dict,
           c("mitosis", "heme", "spermiogenesis", "prolactin",
             "trophoblast", "adipose", "hepatocyte", "follicles",
             "oocyte", "axonemal", "plateletderived", "cuticle",
             "follicle", "centrosome", "morphogenetic", "divide",
             "spermatid", "filopodium", "cameratype", "bile",
             "retina", "cartilage", "notochord", "gut",
             "adipocytes", "myofibril", "pancreas", "fetus",
             "skin", "lamellipodium", "epidermal", "epidermis",
             "kinetochore", "stomach", "drosophila", "muscular",
             "adipocyte", "plus", "keratinocyte", "lactate",
             "muscular", "hematopoietic", "divisions", "flagellar",
             "postembryonic", "digestive", "olfactory", "proliferate",
             "sarcoplasmic", "thrombin", "sprouting", "mesoderm",
             "lymphoid", "hair", "face", "lymph",
             "centromeric", "mucosal", "granulosa", "fibrous",
             "sexual", "pregnancy", "neuroblast", "neuroepithelial",
             "sex", "limb", "male", "uterus",
             "ciliary", "axoneme", "mesenchymal", "flagellum",
             "neuromuscular", "myosin", "chemotaxis", "tooth",
             "marrow", "leukocyte", "taste", "actomyosin",
             "mesenchyme", "ectoderm", "spongiotrophoblast", "placenta",
             "meiotic", "checkpoint", "starvation", "pigment",
             "lungs", "fission", "chromatids", "symbiotic",
             "photoreceptor", "myofibrils", "liver", "platelets",
             "adrenal", "germ", "pharynx", "thyroid",
             "differentiate", "ear", "myeloid", "female",
             "circrna", "radiation", "chromatid", "renal",
             "parathyroid", "mammary", "flavonoid", "keratin",
             "granulocyte", "mesodermal", "sac", "retinoic",
             "melanosomes", "scf", "cytokinesis", "insertion",
             "forebrain", "glands", "centriole", "aorta",
             "intestine", "ganglion", "melanocyte", "wound",
             "dsrna", "differentiated", "smell", "capillary",
             "nephron", "collagen", "retinal", "gland",
             "visual", "somite", "saccharomyces", "angiopoietin",
             "thyrotropinreleasing", "cholangiocytes", "cisretinalretinolbindingprotein", "intraepithelial",
             "virus", "viruses", "hepatoblast", "cisretinol",
             "ears", "diaphragm", "fungi", "hepatocytes",
             "reactioncisretinolretinalbindingprotein", "lectin", "developing", "sclerotome",
             "dermatome", "amalgam", "vertebra", "corneocytes",
             "rankl", "cornification", "embryogenesis", "involucrin",
             "loricrin", "somites", "mycotoxin", "phenylpropanoids",
             "skeleton", "lymphatic", "fungal", "cerevisiae",
             "rectum", "multimerin", "cytolysis", "toxic",
             "protozoa", "dentin", "teeth", "centromere",
             "urothelial", "pelvic", "tailanchored", "globin",
             "urinary", "activin", "peptidoglycan", "hemoglobin",
             "resection", "cornified", "keratinocytes", "vein",
             "mesonephros", "symbionts", "ovulation", "salivary",
             "submandibular", "orthophosphoric", "airway", "bronchi",
             "bronchiole", "ooep", "scmcassociated", "actinmyosin",
             "bacterial", "conifers", "hemicelluloses", "meiosis",
             "myoblast", "muralytic", "transcinnamic", "myotubes",
             "bladder", "rhodopsin", "genitalia", "cilium",
             "fmn", "trachea", "follicular", "tendons",
             "arm", "endocardial", "motor", "mucosa",
             "submucosa", "hindgut", "skeletal", "melanocortin",
             "hemolymph", "gametes", "developmental"),
           c("cyclins", "cyclin", "crosslinking", "crosslink",
             "crosslinks", "transcript", "transcripts", "fibers",
             "fiber"))
           # c("cerevisiae", "saccharomyces", "skeleton", "fiber",
           #   "muscle", "skeletal"),
           # c("endoplasmic"))
           # c("development", "muscle", "differentiation", "blood",
           #   "proliferation", "division", "migration", "cerevisiae",
           #   "morphogenesis"),
           # c("endoplasmic", "extracellular", "atp", "transcription",
           #   "kinase", "synapse", "peptide", "catabolic",
           #   "transmembrane"))
  txt = unlist(strsplit(x, " "))
  txt = Corpus(VectorSource(txt))
  txt = tm_map(txt, PlainTextDocument)
  txt = tm_map(txt, removePunctuation)
  txt = tm_map(txt, removeNumbers)
  txt = tm_map(txt, content_transformer(tolower))
  txt = tm_map(txt, removeWords, c(t.rW))
  # corpus = txt
  # txt = tm_map(txt, stemDocument)
  # txt = tm_map(txt, stemCompletion, corpus)
  tdm = TermDocumentMatrix(txt)
  m = as.matrix(tdm)
  word_freqs = sort(rowSums(m), decreasing=TRUE)
  # word_freqs = word_freqs[word_freqs>1]
  if (is.null(width)) {
    word_freqs = paste(names(word_freqs), collapse=" ")
  } else {
    word_freqs = paste(names(word_freqs)[1:width], collapse="\n")
  }
  gsub("[[:space:]]?NA[[:space:]]?", "", word_freqs)
}
