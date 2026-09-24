# ---------------------------------------------------------------------------
# Separating true allele sequences from PCR/sequencing errors with a graph.
#
# Within one amplicon every distinct read-1 sequence is a vertex; two vertices
# are connected when they differ by exactly one nucleotide (Hamming distance
# 1). Erroneous variants are one substitution away from the true sequence, so
# a true allele sequence is a "hub" with many neighbours ("parent"), while its
# errors are leaves ("children").
# ---------------------------------------------------------------------------

#' Undirected graph connecting sequences at Hamming distance `max_errs`.
sequence_graph <- function(seqs, max_errs = 1) {
  g <- igraph::make_empty_graph(n = length(seqs), directed = FALSE)
  d <- stringdist::stringdistmatrix(seqs, seqs, method = "hamming")
  g <- igraph::add_edges(g, t(which(d == max_errs, arr.ind = TRUE)))
  g <- igraph::simplify(g)
  igraph::set_vertex_attr(g, "label", value = seqs)
}

#' Annotate a cluster table with graph statistics. Adds:
#'   neighbours_read1  number of 1-mismatch neighbours of the read 1
#'   clusters          connected component id
#'   freq              share of the amplicon's reads in this cluster
#'   parents           rank of the vertex inside its component by number of
#'                     neighbours (1 = the most connected sequence)
rank_parents <- function(tab) {
  if (nrow(tab) == 0) return(tab)
  g  <- sequence_graph(tab$read1)
  cc <- igraph::components(g)
  tab$neighbours_read1 <- igraph::degree(g)
  tab$clusters <- cc$membership
  tab$freq     <- prop.table(tab$readnumber)
  parents <- numeric(nrow(tab))
  for (k in unique(cc$membership)) {
    in_k <- tab$clusters == k
    parents[in_k] <- rank(-tab$neighbours_read1[in_k], ties.method = "min")
  }
  tab$parents <- parents
  tab
}
