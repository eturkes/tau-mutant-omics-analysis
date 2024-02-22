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

# This script runs all files in the analysis.
# Sections can be commented out as needed.

setwd(dirname(parent.frame(2)$ofile)) # Move to location of this file.

xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("allen_mouse_hip_ctx_10x.Rmd"),
    output_file = file.path("..", "results", "allen_mouse_hip_ctx_10x.html"),
    envir = new.env()
  )
)

xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("MAPTKI_batch02_01prep.Rmd"),
    output_file = file.path("..", "results", "MAPTKI_batch02_01prep.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("NLGF_MAPTKI_batch02_01prep.Rmd"),
    output_file = file.path("..", "results", "NLGF_MAPTKI_batch02_01prep.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("S305N_batch02_01prep.Rmd"),
    output_file = file.path("..", "results", "S305N_batch02_01prep.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("NLGF_S305N_batch02_01prep.Rmd"),
    output_file = file.path("..", "results", "NLGF_S305N_batch02_01prep.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("P301S_batch02_01prep.Rmd"),
    output_file = file.path("..", "results", "P301S_batch02_01prep.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("NLGF_P301S_batch02_01prep.Rmd"),
    output_file = file.path("..", "results", "NLGF_P301S_batch02_01prep.html"),
    envir = new.env()
  )
)

xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("MAPTKI_batch02_02genefunnel.Rmd"),
    output_file = file.path("..", "results", "MAPTKI_batch02_02genefunnel.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("NLGF_MAPTKI_batch02_02genefunnel.Rmd"),
    output_file = file.path("..", "results", "NLGF_MAPTKI_batch02_02genefunnel.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("S305N_batch02_02genefunnel.Rmd"),
    output_file = file.path("..", "results", "S305N_batch02_02genefunnel.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("NLGF_S305N_batch02_02genefunnel.Rmd"),
    output_file = file.path("..", "results", "NLGF_S305N_batch02_02genefunnel.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("P301S_batch02_02genefunnel.Rmd"),
    output_file = file.path("..", "results", "P301S_batch02_02genefunnel.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("NLGF_P301S_batch02_02genefunnel.Rmd"),
    output_file = file.path("..", "results", "NLGF_P301S_batch02_02genefunnel.html"),
    envir = new.env()
  )
)

xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("comb_batch02_01genefunnel.Rmd"),
    output_file = file.path("..", "results", "comb_batch02_01genefunnel.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("MAPTKI_P301S_batch02_01genefunnel.Rmd"),
    output_file = file.path("..", "results", "MAPTKI_P301S_batch02_01genefunnel.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("NLGF_MAPTKI_NLGF_P301S_batch02_01genefunnel.Rmd"),
    output_file = file.path("..", "results", "NLGF_MAPTKI_NLGF_P301S_batch02_01genefunnel.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("P301S_NLGF_P301S_batch02_01genefunnel.Rmd"),
    output_file = file.path("..", "results", "P301S_NLGF_P301S_batch02_01genefunnel.html"),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("MAPTKI_NLGF_P301S_batch02_01genefunnel.Rmd"),
    output_file = file.path("..", "results", "MAPTKI_NLGF_P301S_batch02_01genefunnel.html"),
    envir = new.env()
  )
)
