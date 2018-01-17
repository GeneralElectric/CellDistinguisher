library(Matrix)                         #Needed by gecd_CellDistinguisher
library(CellMix)                        #Needed by gecd_DeconvolutionCellMix
library(GEOquery)                       #Needed by one-time use of gecd_DataLoader$GSE* routines.
library(biomaRt)

######################################################################
### gecd_CellDistinguisher: finds distinguishers that identify
### specific cell classes without having to know the cell classes.
###
### See the article and supplementary materials for a description of
### the algorithm.
######################################################################

gecd_CellDistinguisher <- function (
    exprLinear,
    genesymb = NULL,
    numCellClasses = 2,
    minDistinguisherAlternatives = 1,
    maxDistinguisherAlternatives = 5,
    minAlternativesLengthsNormalized = 0.7, # 0.0 is no filtering
    expressionQuantileForScale = 0.75,      # 0.0 is no preference for highly expressed genes
    expressionQuantileForFilter = 1.0,      # 1.0 is no filtering.  0.999 is some filtering
    expressionConcentrationRatio = 0.333,   # 0.0 is no filtering
    probesWithGenesOnly = FALSE,
    verbose = 0) {

  ## ##################################################################
  ## The arguments of gecd_CellDistinguisher
  ## ##################################################################
  ##
  ## exprLinear is a matrix of (linear / non-logarithmic) expression
  ## values that has one row for each marker (a.k.a. probe or gene)
  ## and one column for each sample
  ##
  ## numCellClasses is the number of cell classes for which the
  ## function will locate class-specific distinguishers
  ##
  ## minDistinguisherAlternatives is the minimum number of
  ## distinguishers to be located for each cell class
  ##
  ## maxDistinguisherAlternatives is the maximum number of
  ## distinguishers to be located for each cell class
  ##
  ## minAlternativesLengthsNormalized is the threshold for
  ## bestLengthsNormalized values.  Only alternatives with a value at
  ## least as large as the threshold are retained.
  ##
  ## expressionQuantileForScale is in [0, 1].  Low values (e.g., 0.40)
  ## favor distinguishers for a cell class that have very little
  ## expression in the other cell classes.  High values (e.g., 0.80)
  ## shift the balance towards distinguishers with high expression
  ## values.
  ##
  ## genesymb is the associated gene names for the rows of exprLinear.
  ## They are not passed in as names(exprLinear) because multiple rows
  ## (probes) may map to the same gene.
  ##
  ## verbose indicates amount of printed output

  ## ##################################################################
  ## Return values
  ## ##################################################################
  ##
  ## $bestDistinguishers = a matrix of the probes of the discovered
  ## distinguishers.  The first column belongs to the first cell
  ## class, etc.
  ##
  ## $bestLengths = the (diluted) distance of each distinguisher from
  ## the span of the passOneDistinguishers for the other cell classes.
  ##
  ## $bestLengthsNormalized = each $bestLengths value is normalized by
  ## the length for the best distinguisher for its cell subclass
  ##
  ## $passOneDistinguishers = the distinguishers discovered durring
  ## the first pass, which are refined to produce the list of
  ## bestDistinguishers.

  ## ##################################################################
  ## Variables that are global to gecd_CellDistinguisher and its
  ## helper functions
  ## ##################################################################

  ptm <<- proc.time()

  ## ##################################################################
  ## expressionPreprocessing is a helper function for
  ## gecd_CellDistinguisher.
  ## ##################################################################

  expressionPreprocessing <- function (
      exprLinear,
      genesymb,
      expressionQuantileForScale,
      expressionQuantileForFilter,
      expressionConcentrationRatio,
      probesWithGenesOnly,
      verbose) {

    ## ################################################################
    ## Make sure there are no duplicate row names
    ## ################################################################
    if (length(unique(rownames(exprLinear))) != length(rownames(exprLinear))) {
      mesg <- paste("There are only", length(unique(rownames(exprLinear))), "unique row names among the", length(rownames(exprLinear)), "rows of exprLinear; there are duplicates.\n")
      stop(mesg)
    }

    ## ################################################################
    ## If requested, keep only probes that have an associated
    ## genesymb.
    ## ################################################################

    if (probesWithGenesOnly) {
      if (is.null(genesymb)) {
        mesg <- "probesWithGenesOnly=TRUE and genesymb=NULL filters out all probes.\n"
        stop(mesg)
      }
      ProbesWithoutGenes <- (genesymb == "")
      if ((verbose > 3) & (sum(ProbesWithoutGenes) > 0)) {
        cat("Eliminating probes without gene symbols:\n")
        print(rownames(exprLinear)[ProbesWithoutGenes])
      }
      exprLinear <- exprLinear[!ProbesWithoutGenes,]
      genesymb <- genesymb[!ProbesWithoutGenes]
    }

    ## ################################################################
    ## Make sure there are no NA or NaN values
    ## ################################################################

    cNaN <- sum(is.nan(exprLinear))
    cNA <- sum(is.na(exprLinear)) - cNaN # Because each NaN shows up as an NA
    if (cNA + cNaN > 0) {
      mesg <- "The exprLinear matrix includes prohibited"
      if (cNA > 0) {
        mesg <- paste(mesg, "NA")
      }
      if (cNA * cNaN > 0) {
        mesg <- paste(mesg, "and")
      }
      if (cNaN > 0) {
        mesg <- paste(mesg, "NaN")
      }
      mesg <- paste(mesg, "values.\n")
      stop(mesg)
    }

    ## ################################################################
    ## Normalize the sample expressions to put the samples on equal
    ## footing
    ## ################################################################

    ## Note that this normalization has the advantage of effectively
    ## converting FPKM (or RPKM) values to TPM values.

    ## Side note: Arora (2013) assumes that the entries of exprLinear
    ## are integer counts of transcripts (words) and does a minor
    ## correction for the fact that two transcripts selected at random
    ## from a sample are not independent if they are the same instance
    ## of the same transcript.  When the exprLinear values are not
    ## integers (e.g., microarray intensities) this is not applicable.
    ## Even for RNA-seq and integer counts, this correction is de
    ## minimis and is ignored.

    exprLinear <- t(t(exprLinear) / colSums(exprLinear))

    if (verbose > 3) {
      print(proc.time() - ptm); ptm <<- proc.time()
      cat("Columns (samples) normalized\n")
    }

    ## ################################################################
    ## Deal with very high expression values by eliminating probes
    ## with outlier expression levels
    ## ################################################################

    if (expressionQuantileForFilter < 1.0) {
      MyExpressionLimit <- quantile(exprLinear, p=expressionQuantileForFilter)
      ProbesWithExtremeExpression <- apply(exprLinear, 1, function (probe) { return(sum(probe > MyExpressionLimit) > 0) })
      if ((verbose > 3) & (sum(ProbesWithExtremeExpression) > 0)) {
        cat("Eliminating probes with outlier expression:\n")
        print(rownames(exprLinear)[ProbesWithExtremeExpression])
      }
      exprLinear <- exprLinear[!ProbesWithExtremeExpression,]
      if (!is.null(genesymb)) {
        genesymb <- genesymb[!ProbesWithExtremeExpression]
      }
    }

    ## ################################################################
    ## Eliminate those probes where the sample with the highest
    ## expression has much more than the expression level of the
    ## second highest expression value for that probe.
    ## ################################################################

    if (expressionConcentrationRatio > 0.0) {
      ProbesWithConcentratedExpression <- apply(exprLinear, 1, function (probe) { top <- sort(probe, decreasing=TRUE)[1:2] ; return(top[2] <= expressionConcentrationRatio*top[1])})
      if ((verbose > 3) & (sum(ProbesWithConcentratedExpression) > 0)) {
        cat("Eliminating probes with concentrated expression:\n")
        print(rownames(exprLinear)[ProbesWithConcentratedExpression])
      }
      exprLinear <- exprLinear[!ProbesWithConcentratedExpression,]
      if (!is.null(genesymb)) {
        genesymb <- genesymb[!ProbesWithConcentratedExpression]
      }
    }

    ## ################################################################
    ## Make sure the matrix of expression data has enough rank
    ## ################################################################
    matrixRank <- rankMatrix(exprLinear)[1]
    if (matrixRank < numCellClasses) {
      mesg <- sprintf("The expression matrix is of low rank.  Ask for numCellClasses <= %d.\n", matrixRank)
      stop(mesg)
    }
    if (verbose > 3) {
      cat("gecd_CellDistinguisher: rank test complete\n")
    }

    ## ################################################################
    ## Dilute the data
    ## ################################################################

    ## The excess of the expressionQuantileForScale of the exprLinear
    ## values beyond the minimum such value is used as an addend to
    ## the input expression values, to dilute the signal.  Especially
    ## when the quantile is high, the dilution will be more signficant
    ## for low-expression markers, reducing their value as
    ## distinguishers and thus favoring high-expression markers.

    ## This looks for the quantile of the non-zero data.  By focusing
    ## on the non-zero data, we do not get results that depend upon
    ## how many zero-expressed genes were cleaned from the data prior
    ## to this point.  

    aDilution <- quantile(exprLinear[exprLinear > 1e-12], probs = c(0, expressionQuantileForScale), names = FALSE, type= 6)
    dilution <- aDilution[2] - aDilution[1]
    exprLinear <- exprLinear + dilution
    ## Normalize columns (again)
    exprLinear <- t(t(exprLinear) / colSums(exprLinear))

    if (verbose > 2) {
      print(proc.time() - ptm); ptm <<- proc.time()
      cat(sprintf("Expression diluted with Expression[%f quantile] = %g\n", expressionQuantileForScale, dilution))
    }

    ## ################################################################
    ## "Construct" Q, Qbar (pronounced Q-bar), and Qbaradj (Q-bar
    ## adjusted).  The rows of Qbaradj are points in space for which
    ## we want to find the convex hull.
    ##
    ## For speed purposes we will not actually compute
    ## the Q matrices; when we need the Q matrices it is in
    ## combination with other matrices, and the optimal associative
    ## order of multiplying the matrices together in these situations
    ## does not have the computation of Q as an intermediate step.
    ## However, implicitly we will have
    ##           Q = exprLinear       %*% t(exprLinear)
    ##        Qbar = exprLinearBar    %*% t(exprLinear) = Q / rowSums(Q)
    ##     Qbaradj = exprLinearBarAdj %*% t(exprLinear) = t(t(Qbar) - colMeans(Qbar))
    ##
    ## Note that we currently come in with colSums(exprLinear) all
    ## ones, but we will not take advantage of that, just in case this
    ## changes in the future.
    ## ################################################################

    ## We are not explicitly computing Q = exprLinear %*%
    ## t(exprLinear).  We already have exprLinear.
    if (verbose > 3) {
      print(proc.time() - ptm); ptm <<- proc.time()
      cat("Q created, implicitly\n")
    }
    ## Qbar = exprLinearBar %*% t(exprLinear)
    exprLinearBar <- exprLinear / (exprLinear %*% colSums(exprLinear))[, 1]
    if (verbose > 3) {
      print(proc.time() - ptm); ptm <<- proc.time()
      cat("Qbar created, implicitly; i.e., rows of Q normalized\n")
    }
    ## Qbaradj = exprLinearBarAdj %*% t(exprLinear)
    if (TRUE) {
      exprLinearBarAdj <- exprLinearBar   # farthest from origin
    } else {
      exprLinearBarAdj <- t(t(exprLinearBar) - colMeans(exprLinearBar)) # farthest from centroid
      if (verbose > 3) {
        print(proc.time() - ptm); ptm <<- proc.time()
        cat("Qbar centered, implicitly\n")
      }
    }

    ## tExprLinear_exprLinear is a useful intermediate result when the
    ## number of markers for which we have expressions is larger than
    ## the number of samples.  (Otherwise, we should let the
    ## gecd_MatrixChainMultiplication calls that use this product
    ## decide for themselves.)
    tExprLinear_exprLinear <- gecd_MatrixChainMultiplication(t(exprLinear), exprLinear, verbose = verbose - 3)

    return(list(exprLinear = exprLinear,
                tExprLinear_exprLinear = tExprLinear_exprLinear,
                exprLinearBar = exprLinearBar,
                exprLinearBarAdj = exprLinearBarAdj,
                genesymb = genesymb))
  }

  ## ##################################################################
  ## findCandidateDistinguishers is a helper function for
  ## gecd_CellDistinguisher.  It finds an initial set of
  ## distinguishers.
  ## ##################################################################

  findCandidateDistinguishers <- function (exprLinear, tExprLinear_exprLinear, exprLinearBarAdj, numCellClasses) {

    ## ################################################################
    ## selectDistantMarker is a helper function for
    ## findCandidateDistinguishers.  Examining Qbaradj, it finds the
    ## next best distinguisher.
    ## ################################################################

    selectDistantMarker <- function (
        exprLinearBarAdj,
        exprLinear,
        tExprLinear_exprLinear,
        passOneDistinguishers,
        allLengths) {
      ## Compute distances^2 to previous distinguishers at origin
      lengths2 <- gecd_MatrixChainMultiplication(
          exprLinearBarAdj,
          tExprLinear_exprLinear,
          t(exprLinearBarAdj),
          diagonalOnly = TRUE,
          verbose = verbose - 3)
      distinguisher <- which.max(lengths2) # Next distinguisher has maximum distance
      passOneDistinguishers <- c(passOneDistinguishers, distinguisher)
      allLengths <- c(allLengths, sqrt(lengths2[distinguisher]))
      if (verbose > 2) {
        iDistinguisher <- length(passOneDistinguishers)
        print(proc.time() - ptm); ptm <<- proc.time()
        cat(sprintf("First pass: 1 CellClass[%d] distinguisher found = %d \"%s\" at distance %g\n", iDistinguisher, distinguisher, rownames(exprLinear)[distinguisher], allLengths[iDistinguisher]))
      }
      return(list(passOneDistinguishers = passOneDistinguishers, allLengths = allLengths))
    }

    ## ################################################################
    ## Look for the first distinguisher
    ## ################################################################

    passOneDistinguishers <- NULL
    allLengths <- NULL
    if (1 <= numCellClasses) {
      iDistinguisher <- 1

      ans <- selectDistantMarker(
          exprLinearBarAdj,
          exprLinear,
          tExprLinear_exprLinear,
          passOneDistinguishers,
          allLengths)
      passOneDistinguishers <- ans$passOneDistinguishers
      allLengths <- ans$allLengths
      rm(ans)
    }

    ## ################################################################
    ## Look for the second distinguisher
    ## ################################################################

    if (2 <= numCellClasses) {
      iDistinguisher <- 2
      ## Subtract first distinguisher from each distribution
      exprLinearBarAdj <- t(t(exprLinearBarAdj) - exprLinearBarAdj[passOneDistinguishers[iDistinguisher-1], ])

      ans <- selectDistantMarker(
          exprLinearBarAdj,
          exprLinear,
          tExprLinear_exprLinear,
          passOneDistinguishers,
          allLengths)
      passOneDistinguishers <- ans$passOneDistinguishers
      allLengths <- ans$allLengths
      rm(ans)
    }

    ## ################################################################
    ## Look for the third and subsequent distinguishers
    ## ################################################################

    if (3 <= numCellClasses) {
      for (iDistinguisher in 3:numCellClasses) {
        ## Prepare to project new distinguisher to the old
        ## distinguisher (origin)
        project <- exprLinearBarAdj[passOneDistinguishers[iDistinguisher-1], , drop = FALSE]
        ## Useful subexpression
        tExprLinear_exprLinear_tProject <- gecd_MatrixChainMultiplication(
            tExprLinear_exprLinear,
            t(project),
            verbose = verbose - 3)
        projectionNorm <- solve(gecd_MatrixChainMultiplication(
            project,
            tExprLinear_exprLinear_tProject,
            verbose = verbose - 3))
        exprLinearBarAdj <- exprLinearBarAdj - gecd_MatrixChainMultiplication(
            exprLinearBarAdj,
            tExprLinear_exprLinear_tProject,
            projectionNorm,
            project,
            verbose = verbose - 3)
        ans <- selectDistantMarker(
            exprLinearBarAdj,
            exprLinear,
            tExprLinear_exprLinear,
            passOneDistinguishers,
            allLengths)
        passOneDistinguishers <- ans$passOneDistinguishers
        allLengths <- ans$allLengths
        rm(ans)
      }
    }

    return(list(
        passOneDistinguishers = passOneDistinguishers,
        allLengths = allLengths,
        exprLinearBarAdj = exprLinearBarAdj))
  }

  ## ##################################################################
  ## adjustDistinguishers is a helper function for
  ## gecd_CellDistinguisher.  For each distinguisher in the supplied
  ## originalDistinguishers, it projects all the other distinguishers
  ## to the origin and then seeks the list of top
  ## maxDistinguisherAlternatives markers that are most distant from
  ## the origin in the direction of the sole remaining original
  ## distinguisher.  This list likely includes the original
  ## distinguisher, often as the best (most distant) marker.
  ## ##################################################################

  adjustDistinguishers <- function (
      passOneDistinguishers,
      exprLinearBarAdj,
      exprLinear,
      tExprLinear_exprLinear,
      maxDistinguisherAlternatives,
      noDistinguishersYet) {

    ## ################################################################
    ## projectOut is a helper function for adjustDistinguishers.  It
    ## sends multiple distinguishers to the origin.  If there are no
    ## distinguishers yet, then the first distinguisher is sent to the
    ## origin by translation of the space.  All other distinguishers
    ## are sent to the origin by the unique orthogonal projection that
    ## would map them to the origin.
    ## ################################################################

    projectOut <- function (
        passOneDistinguishers,
        exprLinearBarAdj,
        exprLinear,
        tExprLinear_exprLinear,
        noDistinguishersYet) {

      if (length(passOneDistinguishers) > 0 && noDistinguishersYet) {
        ## Subtract first distinguisher from each distribution
        exprLinearBarAdj <- t(t(exprLinearBarAdj) - exprLinearBarAdj[passOneDistinguishers[1], ])
        ## Remove first distinguisher from the vector of original
        ## distinguishers
        passOneDistinguishers <- passOneDistinguishers[-1]
      }
      if (length(passOneDistinguishers) > 0) {
        project <- exprLinearBarAdj[passOneDistinguishers, , drop = FALSE]
        tExprLinear_exprLinear_tProject <- gecd_MatrixChainMultiplication(
            tExprLinear_exprLinear,
            t(project),
            verbose = verbose - 3)
        projectionNormInv <- gecd_MatrixChainMultiplication(
            project,
            tExprLinear_exprLinear_tProject,
            verbose = verbose - 3)
        if (kappa(projectionNormInv) > 1e6) {
          mesg <- sprintf("The 'projectionNormInv' matrix is of low rank (kappa = %g).  Ask for fewer numCellClasses.\n", kappa(projectionNormInv))
          stop(mesg)
        }
        projectionNorm <- solve(projectionNormInv)
        exprLinearBarAdj <- exprLinearBarAdj - gecd_MatrixChainMultiplication(
            exprLinearBarAdj,
            tExprLinear_exprLinear_tProject,
            projectionNorm,
            project,
            verbose = verbose - 3)
      }
      return(list(exprLinearBarAdj = exprLinearBarAdj))
    }

    ## ################################################################
    ## adjust the distinguishers
    ## ################################################################

    stopifnot(nrow(exprLinearBarAdj) == nrow(exprLinear))
    stopifnot(ncol(exprLinearBarAdj) == ncol(exprLinear))
    stopifnot(nrow(tExprLinear_exprLinear) == ncol(exprLinear))
    stopifnot(ncol(tExprLinear_exprLinear) == ncol(exprLinear))

    len <- length(passOneDistinguishers)
    if (len == 1) {
      ## Find the points farthest from the origin in the same
      ## direction as the original distinguisher. 
      project <- exprLinearBarAdj[passOneDistinguishers[1], , drop = FALSE]
      lengths <- gecd_MatrixChainMultiplication(
          exprLinearBarAdj,
          tExprLinear_exprLinear,
          t(project),
          verbose = verbose - 3)
      lengths <- lengths / sqrt(lengths[passOneDistinguishers[1]])
        
      bestDistinguishers <- order(lengths, decreasing = TRUE)[1:maxDistinguisherAlternatives]
      bestLengths <- lengths[bestDistinguishers]
      bestLengthsNormalized <- bestLengths / bestLengths[1]
      if (verbose > 2) {
        print(proc.time() - ptm); ptm <<- proc.time()
        cat(sprintf("Second pass: %d distinguisher(s) found for the cell class with pass-one distinguisher %s (length = %g)\n", maxDistinguisherAlternatives, rownames(exprLinear)[passOneDistinguishers[1]], lengths[passOneDistinguishers[1]]))
        distinguishers <- rownames(exprLinear)[bestDistinguishers]
        if (verbose > 3) {
          print(rbind(distinguishers, bestLengths, bestLengthsNormalized))
        } else {
          print(rbind(distinguishers, bestLengths))
        }
      }
      return(list(
          bestDistinguishers = bestDistinguishers,
          bestLengths = bestLengths,
          bestLengthsNormalized = bestLengthsNormalized))
    } else {

      ## We will process 1:mid1 and mid2:len recursively
      mid1 <- (len + 1) %/% 2
      mid2 <- mid1 + 1
      ## Project out mid2:len then recurse to discovering 1:mid1.
      projectOutExprLinearBarAdj1 <- projectOut(
                                       passOneDistinguishers[mid2:len],
                                       exprLinearBarAdj,
                                       exprLinear,
                                       tExprLinear_exprLinear,
                                       noDistinguishersYet)$exprLinearBarAdj
      ans1 <- Recall(
                passOneDistinguishers[1:mid1],
                projectOutExprLinearBarAdj1,
                exprLinear,
                tExprLinear_exprLinear,
                maxDistinguisherAlternatives,
                noDistinguishersYet = FALSE)
      ## Project out 1:mid1 then recurse to discovering mid2:len.
      projectOutExprLinearBarAdj2 <- projectOut(
                                       passOneDistinguishers[1:mid1],
                                       exprLinearBarAdj,
                                       exprLinear,
                                       tExprLinear_exprLinear,
                                       noDistinguishersYet)$exprLinearBarAdj
      ans2 <- Recall(
                passOneDistinguishers[mid2:len],
                projectOutExprLinearBarAdj2,
                exprLinear,
                tExprLinear_exprLinear,
                maxDistinguisherAlternatives,
                noDistinguishersYet = FALSE)
      return(list(bestDistinguishers = cbind(ans1$bestDistinguishers, ans2$bestDistinguishers),
                  bestLengths = cbind(ans1$bestLengths, ans2$bestLengths),
                  bestLengthsNormalized = cbind(ans1$bestLengthsNormalized, ans2$bestLengthsNormalized)))
    }
  }

  ## ##################################################################
  ## goodEnoughDistinguishers is a helper function for
  ## gecd_CellDistinguisher.  It is responsible for shortening the
  ## list of distinguishers for each cell subtype.
  ## ##################################################################

  goodEnoughDistinguishers <- function (distinguishers, lengths, threshold = 0.5, minLength = 1) {
    ## Remove the distinguishers that do not clear the threshold
    result <- distinguishers
    result[(lengths < threshold) & (row(lengths) > minLength)] <- NA

    ## Remove any duplicate entries except where they are first for a
    ## cell subtype.
    tab <- table(result)
    MyGoodNames <- names(tab)[tab == 1]
    result[!(result %in% MyGoodNames) & !(row(result) == 1)] <- NA
    return(result)
  }

  ## ##################################################################
  ## Do the first pass preprocessing
  ## ##################################################################

  ans <- expressionPreprocessing(exprLinear, genesymb, expressionQuantileForScale, expressionQuantileForFilter, expressionConcentrationRatio, probesWithGenesOnly, verbose)
  exprLinear <- ans$exprLinear
  tExprLinear_exprLinear <- ans$tExprLinear_exprLinear
  exprLinearBar <- ans$exprLinearBar
  exprLinearBarAdj <- ans$exprLinearBarAdj
  genesymb <- ans$genesymb
  rm(ans)

  if (verbose > 1) {
    print(proc.time() - ptm); ptm <<- proc.time()
    cat("gecd_CellDistinguisher: preprocessing complete\n")
  }

  ## ##################################################################
  ## Do the first pass
  ## ##################################################################

  ans <- findCandidateDistinguishers(exprLinear, tExprLinear_exprLinear, exprLinearBarAdj, numCellClasses)
  passOneDistinguishers <- ans$passOneDistinguishers
  allLengths <- ans$allLengths
  rm(ans)

  if (verbose > 1) {
    print(proc.time() - ptm); ptm <<- proc.time()
    cat("gecd_CellDistinguisher: findCandidateDistinguishers complete\n")
  }

  ## ##################################################################
  ## Do the second pass
  ## ##################################################################

  ans <- adjustDistinguishers(
      passOneDistinguishers,
      exprLinearBarAdj,
      exprLinear,
      tExprLinear_exprLinear,
      maxDistinguisherAlternatives,
      noDistinguishersYet = TRUE)
  bestDistinguishers <- ans$bestDistinguishers
  bestLengths <- ans$bestLengths
  bestLengthsNormalized <- ans$bestLengthsNormalized
  rm(ans)

  if (verbose > 1) {
    print(proc.time() - ptm); ptm <<- proc.time()
    cat("gecd_CellDistinguisher: distinguisher adjustment complete\n")
  }

  ## ##################################################################
  ## Do the post processing pass, filtering out poor alternatives
  ## ##################################################################
  bestDistinguishers <- goodEnoughDistinguishers(bestDistinguishers,
                                                 bestLengthsNormalized,
                                                 threshold=minAlternativesLengthsNormalized,
                                                 minLength=minDistinguisherAlternatives)
  bestLengths[is.na(bestDistinguishers)] <- NA
  bestLengthsNormalized[is.na(bestDistinguishers)] <- NA
  ## Bubble the NA values to the lower numbered rows
  if (nrow(bestDistinguishers) > 1) {
    bubbleNA <- function (x) { y <- rep(NA, length(x)) ; z <- x[!is.na(x)] ; if(length(z) > 0) { y[1:length(z)] <- z } ; return(y) }
    bestDistinguishers <- apply(bestDistinguishers, 2, bubbleNA)
    bestLengths <- apply(bestLengths, 2, bubbleNA)
    bestLengthsNormalized <- apply(bestLengthsNormalized, 2, bubbleNA)
  }

  lastGoodRow <- sum(apply(bestDistinguishers, 1, function (rank) {sum(!is.na(rank)) > 0}))
  bestDistinguishers <- bestDistinguishers[1:lastGoodRow, , drop=FALSE]
  bestLengths <- bestLengths[1:lastGoodRow, , drop=FALSE]
  bestLengthsNormalized <- bestLengthsNormalized[1:lastGoodRow, , drop=FALSE]

  bestDistinguishersGeneNames <- NULL
  passOneDistinguishersGeneNames <- NULL
  if (!is.null(genesymb)) {
    passOneDistinguishersGeneNames <- as.character(genesymb[passOneDistinguishers])
    bestDistinguishersGeneNames <- matrix(genesymb[bestDistinguishers], nrow(bestDistinguishers), ncol(bestDistinguishers))
  }
  if (!is.null(rownames(exprLinear))) {
    passOneDistinguishers <- rownames(exprLinear)[passOneDistinguishers]
    bestDistinguishers <- matrix(rownames(exprLinear)[bestDistinguishers], nrow(bestDistinguishers), ncol(bestDistinguishers))
  }

  bestDistinguishersFrame <- data.frame(biolproc=as.vector(col(bestDistinguishers)), rank=as.vector(row(bestDistinguishers)), probeset=as.vector(bestDistinguishers), genesymb=NA, length=as.vector(bestLengths), lengthNormalized=as.vector(bestLengthsNormalized), stringsAsFactors=FALSE)
  passOneFrame <- data.frame(probeset=passOneDistinguishers, genesymb=NA, length=allLengths, stringsAsFactors=FALSE)
  if(!is.null(genesymb)) {
    passOneFrame$genesymb <- passOneDistinguishersGeneNames
    bestDistinguishersFrame$genesymb <- as.vector(bestDistinguishersGeneNames)
  }
  bestDistinguishersFrame <- bestDistinguishersFrame[!is.na(bestDistinguishersFrame$length),]
  rownames(passOneFrame) <- 1:nrow(passOneFrame)
  rownames(bestDistinguishersFrame) <- paste(bestDistinguishersFrame$biolproc, bestDistinguishersFrame$rank, sep="_")

  if (verbose > 1) {
    print(proc.time() - ptm); ptm <<- proc.time()
    cat("gecd_CellDistinguisher: post processing complete\n")
  }

  return(list(
      bestDistinguishers = bestDistinguishers, # deprecated
      bestDistinguishersGeneNames = bestDistinguishersGeneNames, # deprecated
      bestLengths = bestLengths, # deprecated
      bestLengthsNormalized = bestLengthsNormalized, # deprecated
      passOneDistinguishers = passOneDistinguishers, # deprecated
      passOneDistinguishersGeneNames = passOneDistinguishersGeneNames, # deprecated
      passOneLengths = allLengths, # deprecated
      bestDistinguishersFrame = bestDistinguishersFrame,
      passOneFrame = passOneFrame))
}

######################################################################
### gecd_DeconvolutionByDistinguishers: Use the discovered
### distinguishers to deconvolve the original expression data
######################################################################

gecd_DeconvolutionByDistinguishers <- function (
    exprLinear,
    bestDistinguishers,
    nonNegativeOnly = TRUE,
    convexSolution = TRUE,
    verbose = 0) {
  ## exprLinear is the matrix to be factored
  ##
  ## bestDistinguishers is the matrix of anchors that enable the
  ## factorization.  Note that only its first row is used.
  ##
  ## nonNegativeOnly, if set to TRUE strictly enforces that the
  ## entries in the matrix factors be non-negative.
  ##
  ## convexSolution, if set to TRUE strictly enforces that the entries
  ## for a sample in the sampleComposition matrix sum to 1.0.

  ## Even if one or both of nonNegativeOnly and convexSolution are not
  ## TRUE, strong data should produce results reasonably close to the
  ## desired restrictions.

  ## ####################################################################
  ## Pre-process exprLinear parameter.  Make each column sum to 1.
  ## The justifies the convexity requirement for sampleCompositions.
  exprLinear <- t(t(exprLinear) / colSums(exprLinear))

  ## ####################################################################
  ## Pre-process bestDistinguishers parameter: If we were passed all
  ## the bestDistinguishers, not just the best for each cell subclass,
  ## then restrict to the latter for now.
  if (is.matrix(bestDistinguishers)) {
    topDistinguishers <- bestDistinguishers[1,]
  } else {
    topDistinguishers <- bestDistinguishers
  }

  if ((sum(is.na(topDistinguishers)) > 0) | (length(topDistinguishers) != length(unique(topDistinguishers)))) {
    mesg <- "The supplied bestDistinguishers parameter contains NA or duplicates in the first row.\n"
    stop(mesg)
  }

  ## We wish to operate on the undiluted expression data.

  ## QOriginalBar <- exprLinearBar %*% t(exprLinear)
  rowNormalization <- (exprLinear %*% colSums(exprLinear))[, 1]
  exprLinearBar <- exprLinear / rowNormalization
  exprLinearBar[which(rowNormalization == 0.0),] <- 0.0
  tExprLinear_exprLinear <- gecd_MatrixChainMultiplication(
      t(exprLinear),
      exprLinear,
      verbose = verbose - 3)

  ## ####################################################################
  ## Compute cellSubclassSignatures.  
  ## convexCoefficients[probe, cellSubclass] is the
  ## probability of cellSubclass given the probe
  ## ####################################################################

  ## Check the condition of a matrix being inverted
  project <- exprLinearBar[topDistinguishers, , drop=FALSE]
  projectionNormInv <- gecd_MatrixChainMultiplication(
      project,
      tExprLinear_exprLinear,
      t(project),
      verbose = verbose - 3)
  if (kappa(projectionNormInv) > 1e6) {
    mesg <- sprintf("The 'projectionNormInv' matrix is of low rank (kappa = %g).  Ask for fewer numCellClasses.\n", kappa(projectionNormInv))
    stop(mesg)
  }
  projectionNorm <- solve(projectionNormInv)

  ## First, compute convexCoefficientsUnconstrained as those that
  ## approximate the expression data via unconstrained least-squares.
  convexCoefficientsUnconstrained <- gecd_MatrixChainMultiplication(
      exprLinearBar,
      tExprLinear_exprLinear,
      t(project),
      projectionNorm,
      verbose = verbose - 3)

  if (convexSolution) {
    ## Second, adjust the convexCoefficientsUnconstrained so that each
    ## row sums to 1.0.  That is, constrain each probe's approximation
    ## to lie within the hyperplane spanned by the distinguishers.
    constraintCoefficients <- matrix(1.0, 1, length(topDistinguishers))
    constraintConstants <- matrix(1.0, nrow(exprLinear), 1)
    constraintNormInv <- gecd_MatrixChainMultiplication(
        constraintCoefficients,
        projectionNorm,
        t(constraintCoefficients),
        verbose = verbose - 3)
    if (kappa(constraintNormInv) > 1e6) {
      mesg <- sprintf("The 'constraintNormInv' matrix is of low rank (kappa = %g).  Ask for fewer numCellClasses.\n", kappa(constraintNormInv))
      stop(mesg)
    }
    constraintNorm <- solve(constraintNormInv)
    productNorm <- gecd_MatrixChainMultiplication(
        constraintNorm,
        constraintCoefficients,
        projectionNorm,
        verbose = verbose - 3)
    convexCoefficientsConstrained <- convexCoefficientsUnconstrained + (
        gecd_MatrixChainMultiplication(
            constraintConstants,
            productNorm,
            verbose = verbose - 3) -
        gecd_MatrixChainMultiplication(
            convexCoefficientsUnconstrained,
            t(constraintCoefficients),
            productNorm,
            verbose = verbose - 3))
  } else {
    convexCoefficientsConstrained <- convexCoefficientsUnconstrained
  }

  ## Third, enforce non-negativity if requested.
  if (nonNegativeOnly) {
    ## If some coefficients come out negative (by more than epsilon)
    ## then force those to zero and iterate until no remaining
    ## coefficients are forced negative.
    epsilon <- 1e-6
    for (row in which(rowSums(convexCoefficientsConstrained < -epsilon) > 0)) {
      if (convexSolution) {
        ## Start with the "sum to 1" constraint
        constraintCoefficients <- rep(1.0, ncol(convexCoefficientsConstrained))
        constraintConstants <- 1.0
      } else {
        constraintCoefficients <- NULL
        constraintConstants <- NULL
      }
      ## replacement starts out as the value that it will eventually
      ## replace.
      replacement <- convexCoefficientsConstrained[row,]
      while (sum(replacement < -epsilon) > 0) {
        for (col in which(replacement < -epsilon)) {
          ## Add constraints that set a convexCoefficient to zero.
          basis <- rep(0.0, ncol(convexCoefficientsConstrained))
          basis[col] <- 1.0
          constraintCoefficients <- rbind(constraintCoefficients, basis)
          constraintConstants <- cbind(constraintConstants, 0.0)
        }
        ## Compute the constrained convexCoefficients.
        ## It is too slow to check for good conditioning inside this
        ## loop.  Hope for the best!!
        constraintNorm <- solve(gecd_MatrixChainMultiplication(
            constraintCoefficients,
            projectionNorm,
            t(constraintCoefficients),
            verbose = verbose - 3))
        productNorm <- gecd_MatrixChainMultiplication(
            constraintNorm,
            constraintCoefficients,
            projectionNorm,
            verbose = verbose - 3)
        replacement <- convexCoefficientsUnconstrained[row, , drop=FALSE] + (
            gecd_MatrixChainMultiplication(
                constraintConstants,
                productNorm,
                verbose = verbose - 3) -
            gecd_MatrixChainMultiplication(
                convexCoefficientsUnconstrained[row, , drop=FALSE],
                t(constraintCoefficients),
                productNorm,
                verbose = verbose - 3))
      }
      convexCoefficientsConstrained[row,] <- replacement
    }
    ## Everything is now more-or-less positive.  If it is really,
    ## really small, set it to zero.
    convexCoefficientsConstrained[convexCoefficientsUnconstrained < epsilon] <- 0.0
  }

  ## Fourth, apply Bayes' Rule to compute cellSubclassSignatures.
  ## cellSubclassSignatures[probe, cellSubclass] is to be the
  ## probability of probe given the cellSubclass.  We start with
  ## convexCoefficientsConstrained[probe,cellSubclass], which is the
  ## probability of a cellSubclass given a probe.  First multiply its
  ## rows by the prior probability of probes to get joint distribution
  ## over (cellSubclass, probe) pairs.  Second divide by column sums
  ## to get Pr[probe | cellSubclass].

  cellSubclassSignatures <- convexCoefficientsConstrained * (exprLinear %*% colSums(exprLinear))[, 1]
  cellSubclassSignatures <- t(t(cellSubclassSignatures) / colSums(cellSubclassSignatures))

  ## ####################################################################
  ## sampleCompositions[cellSubclass, sample] is the probability of
  ## the cellSubclass given the sample
  ## ####################################################################

  ## Zeroth, restrict cellSubclassSignatures and exprLinear to just
  ## the rows corresponding to distinguishers.  Then renormalize each
  ## column of each to sum to 1.0.

  ## Select a subset of the probes by which to compute sample
  ## compositions.  Are we making the right choice?!!!
  ## goodProbes <- distinguishers[1,]
  goodProbes <- unique(bestDistinguishers[!is.na(bestDistinguishers)])
  ## goodProbes <- rownames(exprLinear)

  exprLinearDistinguishers <- exprLinear[goodProbes, ]
  exprLinearDistinguishers <- t(t(exprLinearDistinguishers) / colSums(exprLinearDistinguishers))
  cellSubclassSignaturesDistinguishers <- cellSubclassSignatures[goodProbes, ]
  cellSubclassSignaturesDistinguishers <- t(t(cellSubclassSignaturesDistinguishers) / colSums(cellSubclassSignaturesDistinguishers))

  ## First, compute sampleCompositionsUnconstrained as those that approximate the
  ## expression data via unconstrained least-squares.
  projectionNormInv <- t(cellSubclassSignaturesDistinguishers) %*% cellSubclassSignaturesDistinguishers
  if (kappa(projectionNormInv) > 1e6) {
    mesg <- sprintf("The 'projectionNormInv' matrix is of low rank (kappa = %g).  Ask for fewer numCellClasses.\n", kappa(projectionNormInv))
    stop(mesg)
  }
  projectionNorm <- solve(projectionNormInv)
  sampleCompositionsUnconstrained <- gecd_MatrixChainMultiplication(
      projectionNorm,
      t(cellSubclassSignaturesDistinguishers),
      exprLinearDistinguishers,
      verbose = verbose - 3)

  if (convexSolution) {
    ## Second, adjust the sampleCompositionsUnconstrained so that each column
    ## (sample's compositions) sums to 1.0.
    constraintCoefficients <- matrix(1.0, length(topDistinguishers), 1)
    constraintConstants <- matrix(1.0, 1, ncol(exprLinear))
    constraintNormInv <- gecd_MatrixChainMultiplication(
        t(constraintCoefficients),
        projectionNorm,
        constraintCoefficients,
        verbose = verbose - 3)
    if (kappa(constraintNormInv) > 1e6) {
      mesg <- sprintf("The 'constraintNormInv' matrix is of low rank (kappa = %g).  Ask for fewer numCellClasses.\n", kappa(constraintNormInv))
      stop(mesg)
    }
    constraintNorm <- solve(constraintNormInv)
    productNorm <- gecd_MatrixChainMultiplication(
        projectionNorm,
        constraintCoefficients,
        constraintNorm,
        verbose = verbose - 3)
    sampleCompositionsConstrained <- sampleCompositionsUnconstrained + (
        gecd_MatrixChainMultiplication(
            productNorm,
            constraintConstants,
            verbose = verbose - 3) -
        gecd_MatrixChainMultiplication(
            productNorm,
            t(constraintCoefficients),
            sampleCompositionsUnconstrained,
            verbose = verbose - 3))
  } else {
    sampleCompositionsConstrained <- sampleCompositionsUnconstrained
  }

  ## Third, enforce non-negativity if requested.
  if (nonNegativeOnly) {
    ## If some composition values come out negative (by more than
    ## epsilon) then force those to zero and iterate until no
    ## remaining coefficients are forced negative.
    epsilon <- 1e-6
    for (col in which(colSums(sampleCompositionsConstrained < -epsilon) > 0)) {
      if (convexSolution) {
        ## Start with the "sum to 1" constraint
        constraintCoefficients <- rep(1.0, nrow(sampleCompositionsConstrained))
        constraintConstants <- 1.0
      } else {
        constraintCoefficients <- NULL
        constraintConstants <- NULL
      }
      ## replacement starts out as the value that it will eventually
      ## replace.
      replacement <- sampleCompositionsConstrained[, col]
      while (sum(replacement < -epsilon) > 0) {
        for (row in which(replacement < -epsilon)) {
          basis <- rep(0.0, nrow(sampleCompositionsConstrained))
          basis[row] <- 1.0
          constraintCoefficients <- cbind(constraintCoefficients, basis)
          constraintConstants <- rbind(constraintConstants, 0.0)
        }
        ## Compute the constrained sampleCompositions.
        ## It is too slow to check for good conditioning inside this
        ## loop.  Hope for the best!!
        constraintNorm <- solve(gecd_MatrixChainMultiplication(
            t(constraintCoefficients),
            projectionNorm,
            constraintCoefficients,
            verbose = verbose - 3))
        productNorm <- gecd_MatrixChainMultiplication(
            projectionNorm,
            constraintCoefficients,
            constraintNorm,
            verbose = verbose - 3)
        replacement <- sampleCompositionsUnconstrained[, col, drop=FALSE] + (
            gecd_MatrixChainMultiplication(
                productNorm,
                constraintConstants,
                verbose = verbose - 3) -
            gecd_MatrixChainMultiplication(
                productNorm,
                t(constraintCoefficients),
                sampleCompositionsUnconstrained[, col, drop=FALSE],
                verbose = verbose - 3))
      }
      sampleCompositionsConstrained[, col] <- replacement
    }
    sampleCompositionsConstrained[sampleCompositionsConstrained < epsilon] <- 0.0
  }

  return(list(
      cellSubclassSignatures = cellSubclassSignatures,
      sampleCompositions = sampleCompositionsConstrained))
}

######################################################################
### gecd_DeconvolutionCellMix is a wrapper for the CellMix
### deconvolution tools.  This is probably most useful with
### method="ssKL" and method="ssFrobenius".
######################################################################

gecd_DeconvolutionCellMix <- function (
                               exprLinear,
                               bestDistinguishers,
                               nonNegativeOnly = TRUE,
                               verbose = 0,
                               log = FALSE,
                               ...) {
  stopifnot(nonNegativeOnly == TRUE)

  MyMarkerList <- NULL
  cCellTypes <- NULL
  if (!is.null(bestDistinguishers)) {
    if (is.vector(bestDistinguishers)) {
      ## Convert vector to a matrix with one row
      bestDistinguishers <- matrix(bestDistinguishers, 1, length(bestDistinguishers))
    }
    ## Convert the matrix to a named list of unamed lists.  There is
    ## one unnamed list per column of the bestDistinguishers matrix.
    MyMarkerList <- split(bestDistinguishers, col(bestDistinguishers))
    ## Name each list by its first distinguisher
    names(MyMarkerList) <- lapply(MyMarkerList, function (x) { return(x[[1]]) } )
    ## Remove any NA entries
    MyMarkerList <- lapply(MyMarkerList, function (x) { return(x[!is.na(x)]) } )
    ## Convert to MarkerList structure
    cCellTypes <- length(MyMarkerList)
    MyMarkerList <- MarkerList(MyMarkerList)
  }

  ## Eliminate probes with no positive expression to avoid confusing
  ## CellMix
  exprLinearReduced <- exprLinear[apply(exprLinear, 1, function (x) { sum(x > 0.0) > 0}),]
  ## Ask CellMix to do its job
  response <- ged(object=ExpressionSet(assayData=exprLinearReduced), x=cCellTypes, data=MyMarkerList, log=log, verbose=verbose, ...)
  ## Reinsert eliminated probes
  cellSubclassSignatures <- matrix(0, nrow(exprLinear), cCellTypes)
  rownames(cellSubclassSignatures) <- rownames(exprLinear)
  colnames(cellSubclassSignatures) <- colnames(response@fit@W)
  cellSubclassSignatures[rownames(response@fit@W),] <- response@fit@W

  return(list(
      cellSubclassSignatures = cellSubclassSignatures,
      sampleCompositions = response@fit@H))
}

######################################################################
### Matrix chain multiplication
######################################################################

gecd_MatrixChainMultiplication <- function (..., verbose = 0, diagonalOnly = FALSE) {

  ## ##################################################################
  ## Variables that are global to this function and its helper
  ## functions
  ## ##################################################################

  ListOfMatrices <- list(...)
  NumMatrices <- length(ListOfMatrices)
  if (verbose > 0) {
    cat(sprintf("Applicable matrix dimensions %s: ", ifelse(diagonalOnly, "(diagonal only)", "")))
    print(c(sapply(ListOfMatrices, nrow), ncol(ListOfMatrices[[NumMatrices]])))
  }
  ## Cost and Best are private member variables that will be used by
  ## the helper routines ComputeOptimum and ComputeProduct.
  ## Cost[i, j] = cost for multiplying matrices i through j, inclusive.
  ## Cost[i, j] is achieved with MatrixProduct[i, k-1] %*% MatrixProduct[k, j]
  Cost <- matrix(0.0, NumMatrices, NumMatrices)
  Best <- matrix(0, NumMatrices, NumMatrices)

  ## ##################################################################
  ## ComputeOptimum is a helper function for
  ## gecd_MatrixChainMultiplication.  It computes the optimal
  ## associative order in which to later multiply the matrices.
  ## ##################################################################

  ComputeOptimum <- function (first, last, diagonalOnly) {
    if (first < last && Best[first, last] == 0) {
      ## We haven't computed this value before.  Let's compute and cache it.
      if (verbose > 1) { cat(sprintf("Computing for range [%d, %d]\n", first, last)) }
      possibilities <- rep(0.0, last-first)
      rowsFirst <- as.numeric(nrow(ListOfMatrices[[first]]))
      if (diagonalOnly) {
        colsLast <- 1
      } else {
        colsLast <- as.numeric(ncol(ListOfMatrices[[last]]))
      }
      for (mid in (first+1):last) {
        rowsMid <- as.numeric(nrow(ListOfMatrices[[mid]]))
        possibilities[mid-first] <-
          Recall(first, mid-1, diagonalOnly = FALSE) +
          Recall(mid, last, diagonalOnly = FALSE) +
          rowsFirst * rowsMid * colsLast
      }
      where <- which.min(possibilities)
      Best[first, last] <<- where + first
      Cost[first, last] <<- possibilities[where]
    }
    ## Return the value from the cache
    return(Cost[first, last])
  }

  ## ##################################################################
  ## ComputeProduct is a helper function for
  ## gecd_MatrixChainMultiplication
  ## ##################################################################

  ComputeProduct <- function (first, last, diagonalOnly) {
    if (first == last) {
      if (verbose > 0) { cat(sprintf("m%d", first)) }
      return(ListOfMatrices[[first]])
    } else {
      mid <- Best[first, last]
      if (verbose > 0) { cat("(") }
      Left <- Recall(first, mid-1, FALSE)
      Right <- Recall(mid, last, FALSE)
      if (diagonalOnly) {
        answer <- rowSums(Left * t(Right))
      } else {
        answer <- Left %*% Right
      }
      if (verbose > 0) { cat(")") }
      return(answer)
    }
  }

  ## ##################################################################
  ## Use ComputeOptimum
  ## ##################################################################

  ComputeOptimum(1, NumMatrices, diagonalOnly)
  if (verbose > 1) {
    print("Best") ; print(Best)
    print("Cost") ; print(Cost)
  }

  ## ##################################################################
  ## Use ComputeProduct
  ## ##################################################################

  answer <- ComputeProduct(1, NumMatrices, diagonalOnly)
  if (verbose > 0) { cat(sprintf(" = %g multiplications\n", Cost[1, NumMatrices])) }
  return(answer)
}

######################################################################
### gecd_mergeListOfDataframes: Use divide and conquer to merge data
### frames relatively quickly.
######################################################################

gecd_mergeListOfDataframes <- function (listOfDataframes) {
  len <- length(listOfDataframes)
  if (len == 1) {
    return(listOfDataframes[[1]])
  } else {
    mid <- len %/% 2
    return(merge(Recall(listOfDataframes[1:mid]), Recall(listOfDataframes[(mid+1):len]), all=TRUE))
  }
}

######################################################################
### gecd_DataLoader$* functions
######################################################################
###
### gecd_DataLoader$* functions return:
###
### exprLinear = linear (not logarithmic) expression values for
### samples from a file or simulation.
###
### genesymb = gene symbols corresponding to probes (rownames) of
### exprLinear, if available
###
### exprLinearClasses = (linear) expression values for cell classes,
### if available

######################################################################
### Structure with no gecd_DataLoader$* functions yet

gecd_DataLoader <- data.frame(nrow=1)[0]   # one row, zero columns

######################################################################
### gecd_DataAnalyzer$* functions:
######################################################################
###
### Attempts at analyzing the data sets
######################################################################

######################################################################
### Structure with no gecd_DataAnalyzer$* functions yet

gecd_DataAnalyzer <- data.frame(nrow=1)[0]   # one row, zero columns

######################################################################
### Enumerate gecd_DataLoader functions

gecd_DataLoader$enumerate.loaders <- function (loaders=gecd_DataLoader, prefix="gecd_DataLoader$") {
  enumerateResponse <- NULL
  for (loader in mixedsort(names(loaders))) {
    if (is.null(names(loaders[[loader]]))) {
      ## This is a leaf node
      if (loader %in% c("enumerate.loaders", "GSETemplate", "help", "mergeListOfDataSets")) {
        ## report nothing
      } else if (loader %in% c("brain.map.H0351", "brain.map.H0351.all", "brain.map.H0351.byregion", "read.table", "GSEGeneric")) {
        ## report the name, but nothing else
        enumerateResponse <- rbind(enumerateResponse, data.frame(loader=paste0(prefix, loader, "()"), numProbes=NA, numSamples=NA, PubMed=NA, minExpression=NA, maxExpression=NA, stringsAsFactors=FALSE))
      } else {
        ## report as much as we can
        z <- loaders[[loader]]()
        enumerateResponse <- rbind(enumerateResponse, data.frame(loader=paste0(prefix, loader, "()"), numProbes=nrow(z$exprLinear), numSamples=ncol(z$exprLinear), PubMed=z$pmid, minExpression=min(z$exprLinear), maxExpression=max(z$exprLinear), stringsAsFactors=FALSE))
      }
    } else {
      ## This is not a leaf node; recurse.
      enumerateResponse <- rbind(enumerateResponse, Recall(loaders[[loader]], paste0(prefix, loader, "$")))
    }
  }
  return(enumerateResponse)
}

gecd_DataLoader$help <- gecd_DataLoader$enumerate.loaders

######################################################################
### Generic read from a file.

gecd_DataLoader$read.table <- function (filename) {
  exprLinear <- data.matrix(read.table(filename, sep="\t", quote="", comment.char="", header=TRUE, row.names=1))
  return(list(
           exprLinear = exprLinear,
           genesymb = NULL,
           exprLinearClasses = NULL,
           pmid = NA))
}


######################################################################
### Simple simulated data

gecd_DataLoader$simulated <- function (numGenes = 22283, numSamples = 21, numCellClasses = 10) {
  exprLog2Classes <- matrix(rnorm(numGenes * numCellClasses, m = 0, sd = 1), ncol = numCellClasses)
  exprLinearClasses <- exprLog2Classes
  exprLinearClasses <- 2^exprLog2Classes
  sampleCompositions <- t(rdirichlet(numSamples, rep(1, numCellClasses)))
  exprLinear <- exprLinearClasses %*% sampleCompositions
  return(list(
           exprLinear = exprLinear,
           genesymb = NULL,
           exprLinearClasses = exprLinearClasses,
           pmid = NA))
}


######################################################################
### A simulated toy example for testing

gecd_DataLoader$toy.example <- function (numCellClasses = 5, numMarkers = 20000, numSamples = 10) {
  ## ####################################################################
  ## concatenateCompositions: create informative sample names
  concatenateCompositions <- function (...) {
    return(paste(substr(paste(round(100+100*..., 0)), 2, 3), collapse="-"))
  }

  ## ####################################################################
  ## makeSymbol: create probe (prefix="p") or gene (prefix="g") names.
  makeSymbol <- function (prefix, index, numDigits) {
    return(paste(prefix, substr(paste(10^numDigits + index), 2, numDigits + 1), sep=""))
  }

  ## ####################################################################
  ## addNoise: The standard deviation of a marker's expression across
  ## the samples is used as a scale for adding noise.  The errorLevel
  ## times this standard deviation is the standard deviation of the
  ## added noise for that marker.

  addNoise <- function (exprLog2, errorLevel = 0.1) {
    xlog.sd <- apply(exprLog2, 1, sd)
    noise <- 0 * exprLog2
    for (i in 1:nrow(exprLog2)) {
      noise[i, ] <- rnorm(ncol(exprLog2), mean = 0, sd = errorLevel*xlog.sd[i])
    }
    return(list(exprLog2 = exprLog2 + noise))
  }

  CellClassExpressions <- matrix(1000*2^rnorm(numMarkers * numCellClasses), numMarkers, numCellClasses)
  ## SampleCompositions <- matrix(rnorm(numCellClasses * numSamples)^2, numCellClasses, numSamples) # Dirichlet(0.5, ...)
  SampleCompositions <- matrix(-log(1-runif(numCellClasses * numSamples)), numCellClasses, numSamples) # Dirichlet(1.0, ...), I think
  ## SampleCompositions[1:numCellClasses, 1:numCellClasses] <- diag(numCellClasses) # identity matrix
  SampleCompositions <- t(t(SampleCompositions) / colSums(SampleCompositions))

  numDigits <- ceiling(log10(numMarkers + 0.5))
  probesymb <- makeSymbol("p", 1:numMarkers, numDigits)
  genesymb <- makeSymbol("g", 1:numMarkers, numDigits)
  numDigits <- ceiling(log10(numCellClasses + 0.5))
  cellClassSymb <- makeSymbol("k", 1:numCellClasses, numDigits)
  sampleSymb <- apply(SampleCompositions, 2, concatenateCompositions)

  rownames(CellClassExpressions) <- probesymb
  colnames(CellClassExpressions) <- cellClassSymb
  rownames(SampleCompositions) <- cellClassSymb
  colnames(SampleCompositions) <- sampleSymb

  SampleExpressions <- CellClassExpressions %*% SampleCompositions

  exprLog2 <- log2(SampleExpressions)
  exprLog2 <- addNoise(exprLog2, 0.1)$exprLog2
  exprLinear <- 2^exprLog2
  return(list(
           exprLinear = exprLinear,
           genesymb = genesymb,
           exprLinearClasses = CellClassExpressions,
           pmid = NA))
}


######################################################################
### GSEGeneric is a function that retrieves and parses a GSE structure
### fetched by getGEO of the GEOquery package.  THESE STRUCTURES ARE
### SOMEWHAT VARIABLE AND THIS FUNCTION DOES NOT ALWAYS WORK.  The
### large size of these data sets can cause the download and parsing
### to be quite slow.  Typical usage is available by running the
### function without arguments.

gecd_DataLoader$GSEGeneric <- function (gseGEO=NULL, retrievedAs=NULL, genesymbColname=NULL, verbose=0) {
  if (is.null(gseGEO) | is.null(retrievedAs)) {
    stop("You must supply the first two arguments.\nTypical usage is:\n  library(GEOquery)\n  options(download.file.method='wget')\n  options(download.file.method.GEOquery='wget')\n  gseGEO <- getGEO('GSE26304', GSEMatrix=FALSE)\n  str(gseGEO@gpls[[1]]@dataTable@table)\n  str(gseGEO@gsms[[1]]@dataTable@columns$Description)\n  MyData <- gecd_DataLoader$GSEGeneric(gseGEO, retrievedAs='log2', genesymbColname='GENE_SYMBOL')\n  MyExprLinear <- MyData[[1]]$exprLinear\n")
  }
  if (verbose >= 1) {
    cat(paste0("Loading data for ", gseGEO@header$geo_accession, "\n"))
  }

  ## Convert the data to a compact format
  platforms <- names(gseGEO@gpls)
  samples <- names(gseGEO@gsms)
  platformsList <- NULL
  for (platform in platforms) {
    if (verbose >= 2) {
      cat(paste0("  Loading data for ", platform, "\n"))
    }
    ## We may have multiple platforms (GPL*).  Process each
    ## separately.
    allSamplesExprLinearDataFrame <- NULL
    allSamplesInfoDataFrame <- NULL
    for (sample in samples) {
      if (gseGEO@gsms[[sample]]@header$platform_id == platform) {
        if (verbose >= 3) {
          cat(paste0("    Loading data for ", sample, "\n"))
        }
        ## This sample is associated with this platform.  Get the
        ## expression values.
        if (verbose >= 4) {
          cat(paste0("      Loading VALUE\n"))
        }
        srcVALUE <- gseGEO@gsms[[sample]]@dataTable@table$VALUE
        if (class(srcVALUE) == "factor") {
          srcVALUE <- levels(srcVALUE)[srcVALUE]
        }
        srcVALUE <- as.numeric(srcVALUE)
        VALUE <- rep(NA, length(srcVALUE))
        if (retrievedAs == "linear") {
          VALUE <- srcVALUE
        } else if (retrievedAs == "log2") {
          VALUE <- 2.0 ^ srcVALUE
        } else if (retrievedAs == "log10") {
          VALUE <- 10.0 ^ srcVALUE
        }
        if (verbose >= 4) {
          cat(paste0("      Loading ID_REF\n"))
        }
        ID_REF <- gseGEO@gsms[[sample]]@dataTable@table$ID_REF
        sampleExprLinearDataFrame <- data.frame(ID_REF, VALUE)
        names(sampleExprLinearDataFrame)[names(sampleExprLinearDataFrame)=="VALUE"] <- sample
        allSamplesExprLinearDataFrame <- c(allSamplesExprLinearDataFrame, list(sampleExprLinearDataFrame))

        if (verbose >= 4) {
          cat(paste0("      Loading sampleInfo\n"))
        }
        ## Load the sample descriptions (sampleInfo)
        sampleInfoDataFrame <- data.frame(SAMPLE=sample)
        sampleInfoRawHeaders <- gseGEO@gsms[[sample]]@header
        for (header in names(sampleInfoRawHeaders)) {
          if (verbose >= 5) {
            cat(paste0("        Loading header=", header, "\n"))
          }
          values <- sampleInfoRawHeaders[[header]]
          for (value in values) {
            split <- strsplit(value, ": *")
            if (length(split[[1]]) == 2) {
              ## This value is a name=value pair
              header <- split[[1]][1]
              value <- split[[1]][2]
            }
            while (header %in% names(sampleInfoDataFrame)) {
              header <- paste0(header, ".")
            }
            sampleInfoDataFrame[1,header] <- type.convert(value, numerals="no.loss")
          }
        }
        rownames(sampleInfoDataFrame)[1] <- sample
        allSamplesInfoDataFrame <- c(allSamplesInfoDataFrame, list(sampleInfoDataFrame))
      }
    }
    ## Now that we have found all samples for this platform, massage
    ## the data structures somewhat
    if (verbose >= 2) {
      cat(paste0("  Massaging exprLinear\n"))
    }
    allSamplesExprLinearDataFrame <- gecd_mergeListOfDataframes(allSamplesExprLinearDataFrame)
    rownames(allSamplesExprLinearDataFrame) <- allSamplesExprLinearDataFrame$ID_REF
    allSamplesExprLinearDataFrame <- allSamplesExprLinearDataFrame[,!(names(allSamplesExprLinearDataFrame) == "ID_REF")]
    exprLinearMatrix <- as.matrix(allSamplesExprLinearDataFrame, drop=FALSE)

    if (verbose >= 2) {
      cat(paste0("  Massaging sampleInfo\n"))
    }
    allSamplesInfoDataFrame <- gecd_mergeListOfDataframes(allSamplesInfoDataFrame)
    rownames(allSamplesInfoDataFrame) <- allSamplesInfoDataFrame$SAMPLE
    allSamplesInfoDataFrame <- allSamplesInfoDataFrame[,!(names(allSamplesInfoDataFrame) == "SAMPLE")]

    ## Load the probe descriptions
    if (verbose >= 2) {
      cat(paste0("  Loading probesetInfo\n"))
    }
    probeNames <- rownames(exprLinearMatrix)
    probesetInfo <- data.frame(ID=probeNames) # rows in same order in probesetInfo and exprLinearMatrix
    probesetInfoRaw <- gseGEO@gpls[[platform]]@dataTable@table
    probesetInfo <- merge(probesetInfo, probesetInfoRaw, all.x=TRUE)
    row.names(probesetInfo) <- probesetInfo$ID

    ## Load the gene symbols (genesymb)
    if (verbose >= 2) {
      cat(paste0("  Loading genesymb\n"))
    }
    genesymb <- NULL
    if (!is.null(genesymbColname)) {
      genesymb <- setNames(probesetInfo[,genesymbColname], row.names(probesetInfo))
    }
    ## Make a onr-element list of a list of useful information for
    ## this platform
    platformList <- list(list(exprLinear=exprLinearMatrix, genesymb=genesymb, sampleInfo=allSamplesInfoDataFrame, probesetInfo=probesetInfo))
    ## Name the one element
    names(platformList) <- list(platform)
    ## Add this platform to the others already computed
    platformsList <- c(platformsList, platformList)
  }
  return(platformsList)
}


######################################################################
### GSE19830
### http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19830
###
### Shen-Orr SS, Tibshirani R, Khatri P, Bodian DL, Staedtler F, Perry
### NM, Hastie T, Sarwal MM, Davis MM, Butte AJ. Cell type-specific
### gene expression differences in complex tissues. Nat Methods. 2010
### Apr;7(4):287-9. doi: 10.1038/nmeth.1439. Epub 2010 Mar 7. PubMed
### PMID: 20208531; PubMed Central PMCID: PMC3699332.
###
### Title: Expression data from pure/mixed brain, liver and lung to
### test feasability and sensitivity of statistical deconvolution
###
### Organism: Rattus norvegicus
###
### Experiment type: Expression profiling by array
###
### Summary: Tissues are often made up of multiple cell-types. Blood,
### for example, contains many different cell-types, each with its own
### functional attributes and molecular signature. In humans, because
### of its accessibility and immune functionality, blood cells have
### been used as a source for RNA-based biomarkers for many
### diseases. Yet, the proportions of any given cell-type in the blood
### can vary markedly, even between normal individuals. This results
### in a significant loss of sensitivity in gene expression studies of
### blood cells and great difficulty in identifying the cellular
### source of any perturbations. Ideally, one would like to perform
### differential expression analysis between patient groups for each
### of the cell-types within a tissue but this is impractical and
### prohibitively expensive.
###
### To test the relationship between measured gene expression in mixed
### samples and the expression of genes in the isolated pure subsets,
### we begin with a situation in which all factors are known. Tissue
### samples from the brain, liver and lung of a single rat were
### analyzed using expression arrays (Affymetrix) in
### triplicate. Homogenates of those three tissues were then mixed
### together at the cRNA level. We then measured the gene expression
### pattern of each mixed sample. Such mixtures mimic the common
### scenario in which biological samples in a dataset are
### heterogeneous and vary in the relative frequency of the component
### subsets from one another.
###
### Overall design: We mixed rat brain, liver and lung biospecimens
### derived from one animal at the cRNA homogenate level in different
### proportions. 3 technical replicates each. Snap frozen rat liver
### and brain was kept frozen while cutting it into pieces. cDNA
### synthesis and labeling was done with a starting amount of 1 μg,
### using the Affymetrix Eukaryotic One-Cycle Target
### Hybridization. Washing, Staining and scanning protocol for
### Eukaryotic Cartridge Arrays with User-Prepared Buffers and
### Sloutions according to the technical manuals (Affymetrix GeneChip
### Expression, Analysis for Cartridge Arrays using the GCAS version
### 1.4, Affymetrix GeneChip Expression Analysis (P/N 701021,Rev. 5) ,
### Affymetrix GeneChip Expression Wash, Stain and Scan (P/N 702731,
### Rev. 3) , following the manufacturer/s instructions. Data was RMA
### normalized.

library(GEOquery)
gecd_DataLoader$GSE19830 <- function () {
  if (FALSE) {
  
    ## Run these steps once in R by cutting and pasting.
    source("CellDistinguisher.R")
    options(download.file.method='wget')
    options(download.file.method.GEOquery='wget')
    GSE19830GEO <- getGEO("GSE19830", GSEMatrix=FALSE)

    ## Examine the output of these two commands for plugging into the third command
    str(GSE19830GEO@gsms[[1]]@dataTable@columns$Description)
    str(GSE19830GEO@gpls[[1]]@dataTable@table)
    GSE19830Data <- gecd_DataLoader$GSEGeneric(GSE19830GEO, retrievedAs="log2", genesymbColname="Gene Symbol")

    ## How many projects did we get?  What do the expression values look like?
    names(GSE19830Data)
    summary(as.vector(GSE19830Data[[1]]$exprLinear))

    ## Keep only first of several alternative gene symbols that are
    ## separated by "//"
    ## GSE19830Data[[1]]$genesymb <- sub("(.*?) //.*", "\\1", GSE19830Data[[1]]$genesymb)

    ## Extract the important parts
    GSE19830ExprLinear <- GSE19830Data[[1]]$exprLinear
    GSE19830genesymb <- GSE19830Data[[1]]$genesymb
    GSE19830SampleInfo <- GSE19830Data[[1]]$sampleInfo
    GSE19830ProbesetInfo <- GSE19830Data[[1]]$probesetInfo

    save(GSE19830ExprLinear, GSE19830genesymb, GSE19830SampleInfo, GSE19830ProbesetInfo, file="GSE19830.RData")
  }

  load("GSE19830.RData")
  return(list(
           exprLinear = GSE19830ExprLinear,
           exprLinearClasses = NULL,
           sampleInfo = GSE19830SampleInfo,
           probesetInfo = GSE19830ProbesetInfo,
           genesymb = GSE19830genesymb,
           pmid = 20208531))
}

gecd_DataAnalyzer$GSE19830 <- function () {
  rm(list=ls())
  source("CellDistinguisher.R")
  MyData <- gecd_DataLoader$GSE19830()
  ## Convert genesymb from R factors to R strings
  MyData$genesymb <- levels(MyData$genesymb)[MyData$genesymb]
  ## Remove probes that are not assoicated with a gene
  MyGoodProbes <- MyData$genesymb != ""
  MyData$exprLinear <- MyData$exprLinear[MyGoodProbes, ]
  MyData$probesetInfo <- MyData$probesetInfo[MyGoodProbes, ]
  MyData$genesymb <- MyData$genesymb[MyGoodProbes]

  if (TRUE) {
    ## Do our best to look at counts of transcripts:
    ## Probes for the same gene will multiply count transcripts for that
    ## gene.  Divide a probe's expression by the number of probes for
    ## its gene.  That way the sum of the probe expression values will
    ## approximate the gene's expression value.
    GenesymbTable <- table(MyData$genesymb)
    MyData$exprLinear <- MyData$exprLinear / as.vector(GenesymbTable[MyData$genesymb])
  }
  ## Compute distinguishers, without using the pure samples
  MyDistinguishers <- gecd_CellDistinguisher(MyData$exprLinear[,-(1:9)], genesymb=MyData$genesymb, numCellClasses=3, probesWithGenesOnly=TRUE)
  ## Change order to liver, brain, lung.  (We know this neworder
  ## because we hand-inspected MyDistinguishers.)
  neworder <- c(1,3,2)
  MyDistinguishers$bestDistinguishers <- MyDistinguishers$bestDistinguishers[, neworder]
  MyDistinguishers$bestDistinguishersGeneNames <- MyDistinguishers$bestDistinguishersGeneNames[, neworder]
  MyDistinguishers$bestLengths <- MyDistinguishers$bestLengths[, neworder]
  MyDistinguishers$bestLengthsNormalized <- MyDistinguishers$bestLengthsNormalized[, neworder]
  MyDistinguishers$passOneDistinguishers <- MyDistinguishers$passOneDistinguishers[neworder]
  MyDistinguishers$passOneDistinguishersGeneNames <- MyDistinguishers$passOneDistinguishersGeneNames[neworder]
  MyDistinguishers$passOneLengths <- MyDistinguishers$passOneLengths[neworder]

  ## Use the distinguishers to deconvolve all pure and mixed samples
  MyDeconvolutionAll <- gecd_DeconvolutionByDistinguishers(MyData$exprLinear, MyDistinguishers$bestDistinguishers, nonNegativeOnly=TRUE, convexSolution=TRUE)
  sampleCompositionsAll <- MyDeconvolutionAll$sampleCompositions
  tissues <- levels(MyData$sampleInfo$tissue)[MyData$sampleInfo$tissue]
  names(tissues) <- MyData$sampleInfo$geo_accession
  colnames(sampleCompositionsAll) <- tissues[colnames(sampleCompositionsAll)]
  rownames(sampleCompositionsAll) <- c("Liver", "Brain", "Lung")
  ## sampleCompositionsAll <- sampleCompositionsAll[, mixedorder(colnames(sampleCompositionsAll))]
  if (FALSE) {
    prmatrix(t(round(sampleCompositionsAll,3)), rowlab=colnames(sampleCompositionsAll), collab=c("liver", "brain", "lung"))
    write.table(t(sampleCompositionsAll), col.names=NA, file="~/Shen-OrrSampleCompositions.txt", sep="\t")
  }

  ## Compare computed sampleCompositionsAll with the sample labels
  liverLabels <- as.numeric(sub("([0-9]+) % Liver.*", "\\1", colnames(sampleCompositionsAll))) / 100.0
  brainLabels <- as.numeric(sub(".*[^0-9]([0-9]+) % Brain.*", "\\1", colnames(sampleCompositionsAll))) / 100.0
  lungLabels <-  as.numeric(sub(".*[^0-9]([0-9]+) % Lung", "\\1", colnames(sampleCompositionsAll))) / 100.0
  sampleCompositionsLabels <- rbind(liverLabels, brainLabels, lungLabels)
  rownames(sampleCompositionsLabels) <- rownames(sampleCompositionsAll)
  colnames(sampleCompositionsLabels) <- colnames(sampleCompositionsAll)

  # Look for the projective transformation that maps the labels to the
  # computed values
  multiplier1 <- median((sampleCompositionsAll[1, 10:42]/sampleCompositionsAll[1, 10:42]) / (sampleCompositionsLabels[1, 10:42]/sampleCompositionsLabels[1, 10:42]))
  multiplier2 <- median((sampleCompositionsAll[2, 10:42]/sampleCompositionsAll[1, 10:42]) / (sampleCompositionsLabels[2, 10:42]/sampleCompositionsLabels[1, 10:42]))
  multiplier3 <- median((sampleCompositionsAll[3, 10:42]/sampleCompositionsAll[1, 10:42]) / (sampleCompositionsLabels[3, 10:42]/sampleCompositionsLabels[1, 10:42]))
  multipliers <- c(multiplier1, multiplier2, multiplier3)
  multipliers <- multipliers * 3 / sum(multipliers)
  tmp1 <- diag(multipliers) %*% sampleCompositionsLabels
  tmp2 <- t(t(tmp1) / colSums(tmp1))    # tmp2 should look a lot like sampleCompositionsAll

  ## Try deconvolution without looking at the pure samples.  Check out
  ## whether it produces reasonable signatures by computing cosine and
  ## correlation.
  MyDeconvolutionMixed <- gecd_DeconvolutionByDistinguishers(MyData$exprLinear[, -(1:9)], MyDistinguishers$bestDistinguishers, nonNegativeOnly=TRUE, convexSolution=TRUE)

  liverIn <- apply(MyData$exprLinear[, 1:3], 1, mean)
  brainIn <- apply(MyData$exprLinear[, 4:6], 1, mean)
  lungIn  <- apply(MyData$exprLinear[, 7:9], 1, mean)
  allIn <- c(liverIn, brainIn, lungIn)

  liverOut <- MyDeconvolutionMixed$cellSubclassSignatures[,1]
  brainOut <- MyDeconvolutionMixed$cellSubclassSignatures[,2]
  lungOut  <- MyDeconvolutionMixed$cellSubclassSignatures[,3]
  allOut <- c(liverOut, brainOut, lungOut)

  cosine <- sum(allIn * allOut) / sqrt(sum(allIn^2) * sum(allOut^2))
  corr <- cor(allIn, allOut)
}


######################################################################
### GSETemplate
### http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSETemplate
###
### INSTRUCTIONS: This assumes a single project (GPL) for the data set.
###
### Change "GSETemplate" to, e.g., "GSE26304"
### Put PMID text citation in header
### Set pmid number in return statment
### Put Title through Overall design sections from GSE page in header
### Set retrievedAs and genesymbColname arguments based on str() commands output

gecd_DataLoader$GSETemplate <- function () {
  if (FALSE) {
    ## Run these steps once in R by cutting and pasting.
    source("../../CellDistinguisher-public.R")
    library("GEOquery")
    options(download.file.method='wget')
    options(download.file.method.GEOquery='wget')
    GSETemplateGEO <- getGEO("GSETemplate", GSEMatrix=FALSE)
    ## Examine the output of these two commands for plugging into the third command
    str(GSETemplateGEO@gsms[[1]]@dataTable@columns$Description)
    str(GSETemplateGEO@gpls[[1]]@dataTable@table)

    GSETemplateData <- gecd_DataLoader$GSEGeneric(GSETemplateGEO, retrievedAs="log2", genesymbColname="GENE_SYMBOL")
    ## How many projects did we get?  What do the expression values look like?
    names(GSETemplateData)
    summary(as.vector(GSETemplateData[[1]]$exprLinear))

    ## Keep only first of several alternative gene symbols that are
    ## separated by "//"
    ## GSETemplateData[[1]]$genesymb <- sub("(.*?) //.*", "\\1", GSETemplateData[[1]]$genesymb)
    ## Extract the important parts
    GSETemplateExprLinear <- GSETemplateData[[1]]$exprLinear
    GSETemplategenesymb <- GSETemplateData[[1]]$genesymb
    GSETemplateSampleInfo <- GSETemplateData[[1]]$sampleInfo
    GSETemplateProbesetInfo <- GSETemplateData[[1]]$probesetInfo
    save(GSETemplateExprLinear, GSETemplategenesymb, GSETemplateSampleInfo, GSETemplateProbesetInfo, file="GSETemplate.RData")
  }

  load("GSETemplate.RData")
  return(list(
           exprLinear = GSETemplateExprLinear,
           exprLinearClasses = NULL,
           sampleInfo = GSETemplateSampleInfo,
           probesetInfo = GSETemplateProbesetInfo,
           genesymb = GSETemplategenesymb,
           pmid = NA))
}

gecd_DataAnalyzer$GSETemplate <- function () {
  source("CellDistinguisher.R")
  GSETemplateData <- gecd_DataLoader$GSETemplate()
}


######################################################################
### gecd_DataLoader$mergeListOfDataSets: Combine a list of data sets
### into one giant data set
######################################################################

gecd_DataLoader$mergeListOfDataSets <- function (listOfDataSets) {
  len <- length(listOfDataSets)
  if (len == 1) {
    return(listOfDataSets[[1]])
  } else {
    mid <- len %/% 2
    Set1 <- Recall(listOfDataSets[1:mid])
    Set2 <- Recall(listOfDataSets[(mid+1):len])

    Set1ExprLinear <- data.frame(Set1$exprLinear)
    Set2ExprLinear <- data.frame(Set2$exprLinear)
    Set1ExprLinear$ID_REF <- row.names(Set1ExprLinear)
    Set2ExprLinear$ID_REF <- row.names(Set2ExprLinear)
    ExprLinear <- merge(Set1ExprLinear, Set2ExprLinear, all=TRUE)
    row.names(ExprLinear) <- ExprLinear$ID_REF
    ExprLinear <- ExprLinear[, !(names(ExprLinear) == "ID_REF")]
    ExprLinear <- ExprLinear[,mixedorder(names(ExprLinear))]
    ExprLinear <- as.matrix(ExprLinear)
    if (length(ExprLinear) == 0) { ExprLinear <- NULL }

    Set1SampleInfo <- data.frame(Set1$sampleInfo)
    Set2SampleInfo <- data.frame(Set2$sampleInfo)
    Set1SampleInfo$SAMPLE <- row.names(Set1SampleInfo)
    Set2SampleInfo$SAMPLE <- row.names(Set2SampleInfo)
    SampleInfo <- merge(Set1SampleInfo, Set2SampleInfo, all=TRUE)
    row.names(SampleInfo) <- SampleInfo$SAMPLE
    SampleInfo <- SampleInfo[, !(names(SampleInfo) == "SAMPLE")]
    SampleInfo <- SampleInfo[mixedorder(row.names(SampleInfo)),]
    if (length(SampleInfo) == 0) { SampleInfo <- NULL }

    Set1ProbesetInfo <- data.frame(Set1$probesetInfo)
    Set2ProbesetInfo <- data.frame(Set2$probesetInfo)
    Set1ProbesetInfo$ID_REF <- row.names(Set1ProbesetInfo)
    Set2ProbesetInfo$ID_REF <- row.names(Set2ProbesetInfo)
    ProbesetInfo <- merge(Set1ProbesetInfo, Set2ProbesetInfo, all=TRUE)
    row.names(ProbesetInfo) <- ProbesetInfo$ID_REF
    ProbesetInfo <- ProbesetInfo[, !(names(ProbesetInfo) == "ID_REF")]
    if (length(ProbesetInfo) == 0) { ProbesetInfo <- NULL }

    Set1Genesymb <- data.frame(GENESYMB = Set1$genesymb)
    Set2Genesymb <- data.frame(GENESYMB = Set2$genesymb)
    Set1Genesymb$ID_REF <- row.names(Set1Genesymb)
    Set2Genesymb$ID_REF <- row.names(Set2Genesymb)
    GenesymbDataFrame <- merge(Set1Genesymb, Set2Genesymb, all=TRUE)
    Genesymb <- as.vector(GenesymbDataFrame$GENESYMB)
    if (length(Genesymb) == 0) {
      Genesymb <- NULL
    } else {
      names(Genesymb) <- GenesymbDataFrame$ID_REF
    }

    return(list(
             exprLinear = ExprLinear,
             exprLinearClasses = NULL,
             sampleInfo = SampleInfo,
             probesetInfo = ProbesetInfo,
             genesymb = Genesymb,
             pmid = NULL))
  }
}

