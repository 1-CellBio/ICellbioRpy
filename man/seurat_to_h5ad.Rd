\name{seurat_to_h5ad}
\alias{seurat_to_h5ad}
\title{Convert Seurat object to h5ad file}
\usage{
seurat_to_h5ad(
  seurat_obj,
  output_file,
  default_assay = "RNA",
  layer = "data",
  include_reductions = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{seurat_obj}{A Seurat object to convert}

\item{output_file}{Character string specifying the output h5ad file path}

\item{default_assay}{Character string specifying which assay to use as the main X matrix (default: "RNA")}

\item{layer}{Character string specifying which layer to use from the assay ("data" or "counts", default: "data")}

\item{include_reductions}{Logical indicating whether to include dimensionality reductions (default: TRUE)}

\item{verbose}{Logical indicating whether to print progress messages (default: TRUE)}
}
\value{
Invisibly returns the output file path
}
\description{
This function converts a Seurat object to an h5ad file format using the anndata package.
It preserves sparse matrix format, metadata, and dimensionality reductions for memory efficiency.
}
\details{
This function converts a Seurat object to h5ad format while preserving:
\itemize{
  \item Sparse matrix format for memory efficiency
  \item Cell metadata from Seurat meta.data
  \item Gene metadata (feature names)
  \item Dimensionality reductions (PCA, UMAP, tSNE, etc.)
  \item Multiple assays as layers in the h5ad file
}

The function uses the anndata package through reticulate to create the h5ad file.
The specified assay and layer will be used as the main X matrix, while other
assays will be stored as additional layers.
}
\examples{
\dontrun{
# Convert Seurat object to h5ad
seurat_to_h5ad(seurat_obj, "output.h5ad")

# Use counts instead of data
seurat_to_h5ad(seurat_obj, "output.h5ad", layer = "counts")

# Specify different default assay
seurat_to_h5ad(seurat_obj, "output.h5ad", default_assay = "SCT")
}
}