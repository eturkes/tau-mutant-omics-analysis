#    This file is part of tau-mutant-omics-analysis.
#    Copyright (C) 2023  Emir Turkes, Naoto Watamura, UK DRI at UCL
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
      if (conserve_memory == TRUE) {
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
      seurat <- RunPCA(seurat, features = rownames(seurat), verbose = FALSE)
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
    mem <- as.numeric(unlist(strsplit(system("free -b", TRUE)[2], " "))[16])
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

word_cloud = function(x, width = NULL){
  t.rW = c("cell", "process", "negative", "positive",
           "activity", "protein", "involved",
           "component", "level", "event", "organismal",
           "cellular", "pathway", "mediated", "dependent",
           "group", "target", "biocarta", "kegg",
           "reactome", "system", "nervous", "cells",
           "time", "structure", "whose", "progression",
           "formation", "divided", "can", "specific",
           "outcome", "two", "form", "one",
           "size", "forms", "becomes", "become",
           "generated", "frequency", "rate", "extent",
           "will", "organism", "inner", "wall",
           "walls", "mammals", "organized", "anatomical",
           "concentration", "directed", "towards", "higher",
           "gradient", "results", "change", "state",
           "terms", "etc", "result", "stops",
           "prevents", "reduces", "activates", "increases",
           "activation", "within", "series", "molecular",
           "modulates", "introducing", "thin", "thick",
           "past", "located", "generation", "molecule",
           "reactions", "pathways", "resulting", "compounds",
           "population", "body", "regulation", "early",
           "commonpartner", "region", "regulated", "means",
           "agent", "structures", "involving", "family",
           "decreases", "together", "set", "daughter",
           "subsequent", "outer", "via", "molecules",
           "comprising", "diameter", "occurring", "removes",
           "another", "adjacent", "presence", "events",
           "specialized", "features", "acquires", "association",
           "type", "include", "consist", "linked",
           "containing", "generate", "member", "widely",
           "usually", "position", "often", "tissue",
           "responsible", "diverse", "range", "leading",
           "contributes", "upper", "part", "slightly",
           "expansion", "internal", "steady", "external",
           "middle", "outflow", "site", "orderly",
           "switching", "also", "known", "effects",
           "defined", "larger", "organisms", "along",
           "powered", "action", "effects", "across",
           "asymmetry", "host", "animal", "pertains",
           "creation", "formed", "tissues", "found",
           "certain", "localized", "distinct", "increase",
           "appearance", "due", "following", "levels",
           "consists", "filters", "end", "intial",
           "either", "double", "introduced", "initiation",
           "addition", "contain", "sorting", "wide",
           "ability", "similar", "act", "primarily",
           "passed", "separation", "received", "steps",
           "initial", "production", "somatic", "shaping",
           "organ", "portion", "systems", "components",
           "key", "mature", "transition", "selection",
           "promotor", "factor", "processing", "cleavage",
           "response", "stimulus", "secretion", "produced",
           "physiologic", "synthesized", "bind", "receptor",
           "trigger", "signaling", "signal", "release",
           "downstream", "communication", "using", "contains",
           "secretion", "transport", "transporter", "movement",
           "response", "signaling", "maintenance",
           "force", "multiplication", "common", "flow",
           "stream", "specialization", "charged", "multicellular",
           "fluid", "proteins", "stem", "converted",
           "transported", "arm", "superfamily", "air",
           "small", "maintained", "parts", "groups",
           "xray", "element", "beam", "strut",
           "rod", "preb", "adaptation", "charges",
           "required", "receive", "characterize", "function",
           "central", "sprouting", "long", "carries",
           "outgoing", "left", "longterm", "neuronal",
           "relatively", "unspecialized", "multiple", "high",
           "number", "tolerance", "induction", "compound",
           "neuron", "work", "perform", "functions",
           "electrical", "biological", "condition", "composed",
           "gives", "neural", "begins", "ends",
           "carry", "innermost", "densely", "packed",
           "mostly", "border", "send", "parallel",
           "brush", "provide", "supplies", "twothirds",
           "principal", "main", "supply", "connecting",
           "commonly", "observed", "visibly", "may",
           "exist", "loosely", "associated", "clusters",
           "combination", "attractive", "fully", "functional",
           "strongly", "several", "enclosed", "individual",
           "give", "rise", "attains", "length",
           "joining", "class", "consisting", "occurs",
           "closely", "related", "may", "many",
           "positions", "classes", "including", "occurs",
           "outside", "location", "connected", "various",
           "decrease", "causes", "makes", "consequence",
           "proceeds", "begins", "circumstances", "require",
           "actual", "numbers", "difficult", "numerous",
           "species", "acts", "upon", "lower",
           "carried", "regular", "longitudinal", "array",
           "assisting", "correct", "pressure", "typically",
           "presentation", "expresses", "passes", "towards",
           "existing", "arising", "present", "side",
           "several", "brown", "polarity", "andor",
           "marginal", "tract", "zone", "indicating",
           "cord", "units", "organic", "chain",
           "distributed", "greater", "without", "concomitant",
           "substance", "dense", "core", "late",
           "substances", "establishment", "initiated", "surface",
           "combining", "involves", "solutes", "role",
           "interaction", "coupled", "classical", "entry",
           "clustering", "domain", "insertion", "disassembly",
           "determination", "sliding", "ion", "potential",
           "acid", "guidance", "assembly", "complex",
           "import", "gland", "fusion", "sequestered",
           "cytosol", "propagation", "channels", "membrane",
           "chemical", "acids", "longchain", "metabolic",
           "unsaturated", "sugar", "carbon", "atoms",
           "base", "signals", "binding", "bonds",
           "derived", "sequestering", "separated", "adaptive",
           "commitment", "signals", "potentials", "extension",
           "ending", "replication", "new", "storage",
           "strands", "sequence", "extrusion", "arrangement",
           "eukaryotic", "organelle", "bounded", "anchored",
           "basally", "elongation", "channel", "smooth",
           "pore", "ions", "systemic", "plasma",
           "anion", "much", "organs", "endogenous",
           "connective", "tuft", "surrounded", "capsule",
           "vertebrate", "intracellular", "precursor", "extracellular",
           "secreted", "cellcell", "proximal", "acute",
           "hydrogen", "alpha", "produces", "transducers",
           "activators", "convey", "cascade", "toward",
           "nuclear", "translocation", "decline", "transduction",
           "accomplished", "retraction", "noncanonical", "bonding",
           "transforming", "right", "ventral", "phase",
           "repulsive", "cues", "canal", "sacs",
           "plant", "natural", "pesticide", "insects",
           "derivative", "linkage", "saturated", "lengths",
           "step", "cycle", "removed", "three",
           "remain", "respectively", "reaction", "ring",
           "transmission", "increased", "white", "achieved",
           "decreased", "attachment", "anterior", "secretory",
           "temporary", "negatively", "nucleus", "droplets",
           "residues", "linkages", "bonding", "incorporation",
           "context", "transfer", "multisubunit", "subunits",
           "processes", "identical", "sister", "environment",
           "density", "detection", "stimuli", "exposure",
           "soluble", "solvents", "powers", "defense",
           "immature", "break", "coat", "exogenous",
           "origin", "tertiary", "cycles", "novo",
           "always", "dorsal", "processes", "nuclear",
           "transmission", "loss", "corpus", "removing",
           "input", "acquire", "modulation", "generally",
           "hemispheres", "covalently", "bonded", "bonding",
           "complete", "apical", "protrusion", "subcellular",
           "compartment", "striated", "examples", "driven",
           "doublet", "tube", "planar", "closure",
           "coordinated", "plane", "partly", "axis",
           "activated", "removal", "includes", "identity",
           "sheath", "stage", "occur", "chains", "cushion",
           "glands", "focal", "stages", "evoked",
           "foam", "confining", "lateral", "cavity",
           "third", "interactions", "node", "mediates",
           "surroundings", "cyclic", "intermediate", "mass",
           "selfpropelled", "bound", "excluding", "point",
           "modification", "particle", "reactive", "active",
           "segregation", "substrate", "gamma", "execution",
           "facial", "branches", "fibers", "metal",
           "residue", "radial", "satellite", "restore",
           "interacting", "selectively", "noncovalently", "biosynthetic",
           "inheritance", "links", "light", "link",
           "tight", "enables", "innate", "biosynthesis",
           "tubule", "bulb", "cohesion", "replicated",
           "lowdensity", "directly", "solute", "lens",
           "human", "roof", "poles", "structural",
           "gene", "specification", "allows", "direct",
           "exchange", "smoothened", "secondary", "controlled",
           "intrinsic", "relaxation", "strength", "animals",
           "bodies", "information", "bond", "products",
           "embedded", "direction", "characteristics", "rhythm",
           "changes", "major", "requires", "residing",
           "physiological", "capable", "mediator", "amounts",
           "basement", "leftright", "pattern", "patterns",
           "twelve", "general", "transmitting", "plays",
           "main", "maintains", "providing", "pathwayrestricted",
           "targeting", "creating", "acquisition", "free",
           "linear", "branching", "interpreted", "contribute",
           "unit", "helps", "repeating", "constituent",
           "progresses", "thus", "surfaces", "produce",
           "area", "develop", "covers", "consequent",
           "conditions", "mixture", "affects", "families",
           "normal", "favoring", "replacement", "around",
           "rid", "play", "retains", "throughout",
           "life", "highly", "ordered", "generates",
           "combinations", "adapts", "contributing", "emerge",
           "vessel", "preexisting", "travel", "lies",
           "called", "example", "four", "space",
           "first", "initiate", "characteristic", "average",
           "predominantly", "alternating", "second", "opens",
           "actions", "single", "forming", "directing",
           "contained", "matter", "naturally", "possessing",
           "center", "help", "exhibit", "preinitiation",
           "machinery", "abortive", "repeatedly", "attached",
           "receives", "separates", "liquid", "material",
           "independent", "tail", "latent", "barrier",
           "impulse", "distal", "substituted", "least",
           "according", "neutral", "travels", "brain",
           "ear", "partitioning", "choice", "source",
           "III", "blos", "coloring", "insecticide",
           "volumesensitive", "box", "pigment", "bundle",
           "chctype", "outward", "divergent", "conserved",
           "driving", "finger", "organization", "spontaneous",
           "envelope", "depolarization", "minusend", "comprises",
           "proteincoupled", "ruffle", "spindle", "cortex",
           "thiolester", "rnan", "rcosr", "localization",
           "dorsalventral", "segment", "granules", "sound",
           "acceptor", "primary", "network", "mechanism",
           "bubble", "vibration", "trimers", "domains",
           "activating", "heavy", "different", "retention",
           "apparatus", "pad", "tubular", "organizing",
           "faces", "subunit", "lungs", "sshaped",
           "binds", "nonhomologous", "steadystate", "hap",
           "responsive", "latency", "trunk", "segments",
           "iris", "adaptor", "innervates", "touch",
           "basal", "egg", "contents", "vitamin",
           "anions", "thalamus", "titin", "regulator",
           "posture", "controlling", "voluntarily", "alignment",
           "placode", "half", "monoatomic", "cardioblast",
           "lymphoid", "morphogenesis", "essential", "mechanisms",
           "membranebounded", "absence", "atom", "genes",
           "cerevisiae", "axons", "acting", "absorption",
           "budding", "centromeric", "acetate", "appears",
           "receptors", "coenzyme", "additional", "acidresponsive",
           "dinucleotide", "beta", "areas", "alphabeta",
           "activator", "characterized", "approximately", "covalent",
           "map", "basic", "axonal", "catalyzes",
           "repeat", "repeats", "articular", "electrondense",
           "enzyme", "consensus", "associate", "minusenddirected",
           "extrinsic", "branch", "klinked", "bud")
  txt = unlist(strsplit(x, " "))
  txt = Corpus(VectorSource(txt))
  txt = tm_map(txt, PlainTextDocument)
  txt = tm_map(txt, removePunctuation)
  txt = tm_map(txt, removeNumbers)
  txt = tm_map(txt, content_transformer(tolower))
  txt = tm_map(txt, removeWords, c(t.rW, stopwords("english"), stopwords("SMART")))
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
