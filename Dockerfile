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

FROM rocker/rstudio:4.3.2

LABEL org.opencontainers.image.authors="Emir Turkes emir.turkes@eturkes.com"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        libglpk-dev \
        libudunits2-0 \
        libproj22 \
        libgdal30 \
        liblzma-dev \
        libbz2-dev \
    && Rscript -e "install.packages('rmarkdown')" \
        -e "install.packages('markdown')" \
        -e "install.packages('conflicted')" \
        -e "install.packages('DT')" \
        -e "install.packages('Seurat')" \
        -e "install.packages('viridis')" \
        -e "install.packages('factoextra')" \
        -e "install.packages('pals')" \
        -e "install.packages('tm')" \
        -e "install.packages('qgraph')" \
        -e "install.packages('fastmatch')" \
        -e "install.packages('BiocManager')" \
        -e "install.packages('remotes')" \
        -e "BiocManager::install('SingleCellExperiment')" \
        -e "BiocManager::install('scuttle')" \
        -e "BiocManager::install('glmGamPoi')" \
        -e "BiocManager::install('biomaRt')" \
        -e "BiocManager::install('limma')" \
        -e "BiocManager::install('GSEABase')" \
        -e "BiocManager::install('GSVA')" \
        -e "BiocManager::install('edgeR')" \
        -e "BiocManager::install('ComplexHeatmap')" \
        -e "BiocManager::install('GO.db')" \
        -e "BiocManager::install('DropletUtils')" \
        -e "BiocManager::install('scDblFinder')" \
        -e "remotes::install_github('jmw86069/multienrichjam', dependencies = TRUE)" \
        -e "remotes::install_github('immunogenomics/presto')" \
    && apt-get clean \
    && rm -Rf /var/lib/apt/lists/ \
        /tmp/downloaded_packages/ \
        /tmp/*.rds
