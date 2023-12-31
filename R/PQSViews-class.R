###
## PQSViews class
##
## Author: Jiri Hon <jiri.hon@gmail.com>
## Date: 2016/01/17
## Package: pqsfinder
##


#' An S4 class to represent potential quadruplex forming sequences
#'
#' Represents potential quadruplex forming sequences found by
#' \code{\link{pqsfinder}} function. This is a subclass of
#' \code{\link{XStringViews-class}} class and adds one more slot.
#'
#' @slot density Numbers of PQS (potential quadruplex forming sequences)
#'               overlapping at each position in input sequence.
#' @slot max_scores Score of the best PQS found at each position.
#'
.PQSViews <- setClass(
  "PQSViews",
  contains = "XStringViews",
  slots = c(
    density = "integer",
    max_scores = "integer"
  ),
  validity = function(object) {
    if (!is.na(object@density[1]) && length(object@subject) != length(object@density) ) {
      return("Length of the density vector is not equal to the length of the subject.")
    }
    if (!is.na(object@max_scores[1]) && length(object@subject) != length(object@max_scores)) {
      return("Length of the score distribution vector is not equal to the length of the subject.")
    }
    return(TRUE)
  }
)


#' PQSViews class constructor
#'
#' User friendly constructor for PQSViews class representing potential
#' quadruplex forming sequences (PQS). PQSViews is a subclass of
#' \code{\link{XStringViews}} class and adds two more slots to store
#' PQS density and PQS score distribution.
#' 
#' Use \code{\link{elementMetadata}} function to get extra PQS features
#' like number of tetrads, bulges, mismatches or loop lengths.
#'
#' @param subject DNAString or RNAString object.
#' @param start Start positions.
#' @param width Lengths.
#' @param strand Strand specifications.
#' @param score Scores.
#' @param density Numbers of PQS overlapping at each position in \code{subject}.
#' @param max_scores Score of the best PQS found at each position.
#' @param nt Tetrad numbers.
#' @param nb Bulge counts.
#' @param nm Mismatch counts.
#' @param rl1 Run 1 lengths.
#' @param rl2 Run 2 lengths.
#' @param rl3 Run 3 lengths.
#' @param ll1 Loop 1 lengths.
#' @param ll2 Loop 2 lengths.
#' @param ll3 Loop 3 lengths.
#' @return PQSViews object.
#'
#' @examples
#' pv <- PQSViews(DNAString("GGTGGTGGTGG"), 1, 11, "+", 33, as.integer(rep(1, 11)),
#'                as.integer(rep(33, 11)), 2, 0, 0, 2, 2, 2, 1, 1, 1)
#' start(pv)
#' width(pv)
#' strand(pv)
#' score(pv)
#' density(pv)
#' maxScores(pv)
#' elementMetadata(pv)
#'
PQSViews <- function(
  subject, start, width, strand, score, density, max_scores,
  nt, nb, nm, rl1, rl2, rl3, ll1, ll2, ll3)
{
  ix <- order(start)
  .PQSViews(
    subject = subject, ranges = IRanges(start = start[ix], width = width[ix]),
    elementMetadata = DataFrame(
      strand = strand[ix], score = score[ix], nt = nt[ix],
      nb = nb[ix], nm = nm[ix],
      rl1 = rl1[ix], rl2 = rl2[ix], rl3 = rl3[ix],
      ll1 = ll1[ix], ll2 = ll2[ix], ll3 = ll3[ix]
    ),
    density = density, max_scores = max_scores
  )
}


#' Get PQS score vector
#'
#' @param x PQSViews object.
#' @return Score vector.
#' @examples
#' pqs <- pqsfinder(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
#' score(pqs)
#'
setMethod("score", "PQSViews", function(x) mcols(x)$score)

#' Get PQS strand vector
#'
#' @param x PQSViews object.
#' @return Strand vector.
#' @examples
#' pqs <- pqsfinder(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
#' strand(pqs)
#'
setMethod("strand", "PQSViews", function(x) mcols(x)$strand)

#' Get density vector
#'
#' Desity vector represents numbers of PQS (potential quadruplex forming
#' sequences) overlapping at each position in input sequence.
#'
#' @param x PQSViews object.
#' @return Density vector.
#' @examples
#' pqs <- pqsfinder(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
#' density(pqs)
#'
setMethod("density", "PQSViews", function(x) x@density)

#' Get vector of maximal scores
#' 
#' Get vector of maximal scores for a given object.
#' 
#' @param x An object.
#' @param ... Additional arguments, for use in specific methods.
#' @return Vector of maximal scores.
#' @examples
#' showMethods("maxScores")
#' 
setGeneric("maxScores", function(x, ...) {
  standardGeneric("maxScores")
})

#' Get vector of maximal scores
#'
#' For each sequence position it gives the maximal score
#' of all PQS conformations which overlap that position.
#' 
#' @param x PQSViews object.
#' @return Vector of maximal scores.
#' @examples
#' pqs <- pqsfinder(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
#' maxScores(pqs)
#'
setMethod('maxScores', 'PQSViews', function(x) x@max_scores)

## The 2 helper functions below convert a given view on an XString object
## into a character-string.
##
## Both assume that 'start' <= 'end' (so they don't check it) and
## padd the result with spaces to produce the "margin effect"
## if 'start' or 'end' are out of limits.
##
## NOTE: Heavily based on Biostrings package, file XStringViews-class.R
##

## nchar(get_view(x, start, end)) is always end-start+1
##
.get_view <- function(x, start, end)
{
  as.character(subseq(x, start, end))
}

## nchar(get_snippet(x, start, end, snippetWidth)) is <= snippetWidth
##
.get_snippet <- function(x, start, end, snippetWidth, strand)
{
  if (snippetWidth < 7)
    snippetWidth <- 7
  width <- end - start + 1
  if (width <= snippetWidth) {
    .get_view(x, start, end)
  } else {
    w1 <- (snippetWidth - 2) %/% 2
    w2 <- (snippetWidth - 3) %/% 2
    paste(.get_view(x, start, start+w1-1),
          "...",
          .get_view(x, end-w2+1, end), sep="")
  }
}


## Show header of output table
##
.show_vframe_header <- function(iW, cols)
{
  cat(format("", width=iW)) # Print padding

  for (col in cols)
  {# Print column names
    cat(" ")
    cat(format(col$nm, width=col$width, justify="right"))
  }
  cat("\n")
}


## Show row of output table
##
.show_vframe_line <- function(x, i, iW, cols)
{
  # Print PQS index
  cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"))

  colW <- 0 # Sum of all column width

  for (col in cols)
  {# Print column values
    cat(" ")
    value <- do.call(col$fn, list(x))[i]
    cat(do.call("format", c(value, col[3:length(col)], justify="right")))

    colW = colW + col$width
  }
  snippetW <- getOption("width") - iW - colW - length(cols) - 3
  cat(" [",
      .get_snippet(subject(x), start(x)[i], end(x)[i], snippetW),
      "]\n", sep="")
}


## Shot dots in place of hidden table rows
##
.show_vframe_line_dots <- function(x, iW, cols)
{
  cat(format("...", width=iW, justify="right"))

  for (col in cols) {
    cat(" ")
    cat(format("...", width=col$width, justify="right"))
  }
  cat("\n")
}

.get_nt <- function(pv) elementMetadata(pv)$nt
.get_nb <- function(pv) elementMetadata(pv)$nb
.get_nm <- function(pv) elementMetadata(pv)$nm
.get_rl1 <- function(pv) elementMetadata(pv)$rl1
.get_rl2 <- function(pv) elementMetadata(pv)$rl2
.get_rl3 <- function(pv) elementMetadata(pv)$rl3
.get_ll1 <- function(pv) elementMetadata(pv)$ll1
.get_ll2 <- function(pv) elementMetadata(pv)$ll2
.get_ll3 <- function(pv) elementMetadata(pv)$ll3

## Show all output table rows
## 'half_nrow' must be >= 1
##
.show_vframe <- function(x, half_nrow=9L)
{
  ## Column definitions
  ## nm = Column header
  ## fn = Function name to get row values
  ## ... other parameters passed to format function
  ##
  cols <- list(
    list(nm="start",  fn="start" ),
    list(nm="width",  fn="width" ),
    list(nm="score",  fn="score" ),
    list(nm="strand", fn="strand"),
    list(nm="nt", fn=".get_nt"),
    list(nm="nb", fn=".get_nb"),
    list(nm="nm", fn=".get_nm")
  )

  i <- 1
  for (col in cols)
  {# Calculate column widths
    col_max <- max(do.call(col$fn, list(x)))
    col_maxstr <- do.call("format", c(col_max, col[3:length(col)]))
    cols[[i]] <- c(col, width=max(nchar(col_maxstr), nchar(col$nm)))
    i = i + 1
  }

  lx <- length(x)
  iW <- nchar(format(lx)) + 2 # Two extra for square brackets

  .show_vframe_header(iW, cols)

  if (lx <= 2*half_nrow + 1)
  {# Show all
    for (i in seq_len(lx))
      .show_vframe_line(x, i, iW, cols)
  }
  else
  {# Show first and last views
    for (i in 1:half_nrow)
      .show_vframe_line(x, i, iW, cols)

    .show_vframe_line_dots(x, iW, cols)

    for (i in (lx-half_nrow+1L):lx)
      .show_vframe_line(x, i, iW, cols)
  }
}


#' Show method
#'
#' @param object PQSViews object.
#' @return PQSViews object printed.
#'
setMethod("show", "PQSViews", function(object)
{
  subject <- subject(object)
  lsub <- length(subject)

  cat("  PQS views on a ", lsub, "-letter ", class(subject),
      " subject", sep="")
  cat("\nsubject:", .get_snippet(subject, 1, lsub, getOption("width")-9, "+"))
  cat("\nquadruplexes:")

  if (length(object) == 0) {
    cat(" NONE\n")
  }
  else {
    cat("\n")
    .show_vframe(object)
  }
})


## Set names for XStringSet
##
.get_xstring_names <- function(from) {
  xstring_names <- sprintf(
    "pqsfinder;G_quartet;start=%s;end=%s;strand=%s;score=%s;nt=%s;nb=%s;nm=%s;rl1=%s;rl2=%s;rl3=%s;ll1=%s;ll2=%s;ll3=%s;",
    start(from), end(from), strand(from), score(from),
    .get_nt(from), .get_nb(from), .get_nm(from),
    .get_rl1(from), .get_rl2(from), .get_rl3(from),
    .get_ll1(from), .get_ll2(from), .get_ll3(from)
  )
  return(xstring_names)
}


###
## Coerce PQSViews to DNAStringSet
##
setAs("PQSViews", "DNAStringSet", function(from)
{
  s <- DNAStringSet(subject(from), start(from), end(from))
  names(s) <- .get_xstring_names(from)
  return(s)
})


###
## Coerce PQSViews to RNAStringSet
##
setAs("PQSViews", "RNAStringSet", function(from)
{
  s <- RNAStringSet(subject(from), start(from), end(from))
  names(s) <- .get_xstring_names(from)
  return(s)
})


###
## Coerce PQSViews to GRanges
##
setAs("PQSViews", "GRanges", function(from)
{
  source <- rep("pqsfinder", length(from))
  type <- rep("G_quartet", length(from))
  seqnames <- rep("chr1", length(from))
  seqlen <- length(subject(from))
  names(seqlen) <- "chr1"
  
  GRanges(
    seqnames,
    IRanges(start(from), end(from)),
    strand(from),
    score = score(from),
    seqlengths = seqlen,
    nt = .get_nt(from),
    nb = .get_nb(from),
    nm = .get_nm(from),
    rl1 = .get_rl1(from),
    rl2 = .get_rl2(from),
    rl3 = .get_rl3(from),
    ll1 = .get_ll1(from),
    ll2 = .get_ll2(from),
    ll3 = .get_ll3(from),
    source = source,
    type = type
  )
})


#' Coerce to character vector
#'
#' @param x PQSViews object.
#' @return Character vector representing PQS.
#'
setMethod("as.character", "PQSViews", function(x)
{
  if (class(subject(x)) == "DNAString") {
    s <- as(x, "DNAStringSet")
  } else {
    s <- as(x, "RNAStringSet")
  }
  as.character(s)
})


###
## Convert to printable string
##
setMethod("toString", "PQSViews", function(x, ...)
{
  toString(as.character(x), ...)
})
