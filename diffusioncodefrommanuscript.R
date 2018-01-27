
diffuseMat <- function(data, ndims = 4, nn = 0.2, sigma = 12,
                       removeFirst = TRUE, bluntDensityCorrection = FALSE,
                       useARPACK = TRUE, normBy = c("Column", "Row")) {
  normBy <- tolower(match.arg(normBy))
  nnp <- nn
  nn <- ceiling(ncol(data) * nn)      # Number of nearest neighbours to include
  KT <- sigma^2                       # Diffusion scale parameter
  cat("Calculating distance matrix. Please wait...")
  d2 <- as.matrix(fastDist(data, squared = TRUE))       # distance matrix calculation
  cat("Done.\n")
  cat("Applying Gaussian kernel.")
  W <- exp(-d2 / (2*KT))              # Apply Gaussian kernel
  if (nnp < 1 & nnp > 0) {
    cat("Calculating and retaining nearest neighbours.\n")
    R <- apply(d2, 2, function(x) sort(x)[nn])
    R <- matrix(rep(R,ncol(d2)), ncol = ncol(d2)) # Find distance for nn closest neighbour
    W <- (d2<R) * W                     # Only keep nn nearest neighbours
  }
  cat("Calculating and applying local density correction.\n")
  D <- colSums(W)
  q <- D %*% t(D)                     # Calculate local density
  if (bluntDensityCorrection) q <- sqrt(q)
  diag(W) <- 0
  H <- W / q                          # Correct for local density
  colS <- colSums(H)
  rowS <- rowSums(H)
  if (normBy[1]=="column") {
    eS <- colS == 0                 # Excluded cells with no transitions
    Hp <- t(t(H[!eS,!eS]) / colS[!eS])  # Normalise matrix
  } else if (normBy[1]=="row") {
    eS <- rowS == 0
    Hp <- H[!eS,!eS] / rowS[!eS]
  } else if (normBy[1]=="original") {
    eS <- colS == 0
    Hp <- H[!eS,!eS] / colS[!eS]
  } else cat("No normalistion performed.\n")
  cat("Calculating eigen vectors. Please wait...")
  n <- nrow(d2)
  chooseDims <- function(E) {     # Sub-function to sort and select largest eigenvecs/vals
    if (ncol(E$vectors) <= ndims) ndims <- ncol(E$vectors) - 1
    eigOrd <- order(Re(E$values), decreasing = TRUE)
    E$values <- Re(E$values[eigOrd][startDim:(ndims+1)])
    E$vectors <- Re(E$vectors[,eigOrd][,startDim:(ndims + 1)]) # Remove first eigen vector
    rownames(E$vectors) <- colnames(data)[!eS]
    colnames(E$vectors) <- 1:ncol(E$vectors)
    return(E)
  }
  if (useARPACK) {
    decomp <- eigs(Hp, which = "LR", k = ndims + 1)
    if (removeFirst) {
      startDim <- 2
    } else startDim <- 1
    decomp <- chooseDims(decomp)
    decomp$usedARPACK <- TRUE
    if (length(which(eS)) != 0) cat(paste("\nWarning:\nCells \"", paste(names(which(eS)), collapse = ", "), "\" have no transitions and have been removed from the analysis\n", sep=""))
  } else {
    decomp <- eigen(Hp)                      # Eigen decomposition
    if (removeFirst) {
      startDim <- 2
    } else startDim <- 1
    decomp <- chooseDims(decomp)
    decomp$nconv <- integer(0)
    decomp$niter <- integer(0)
    decomp$usedARPACK <- FALSE
  }
  decomp$nn <- nnp
  return(decomp)
}

fastDist <- function(x, squared = FALSE) {
  a <- colSums(x^2)
  a <- a * matrix(1, ncol(x), ncol(x))
  a <- a + t(a)
  ab <- t(x) %*% x
  d <- a - 2 * ab
  diag(d) <- 0
  if (!squared) d <- sqrt(d)
  return(d)
}

