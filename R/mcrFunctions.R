# Functions that support the cghMCR class
#
# Copyright 2006, J. Zhang, all rights reserved
#

# Gets the MCRs based on the segmentation data contained by segList.
# segList is a data frame extracted from the output of DNAcopy
getMCR <- function(cghmcr){
  
  getMCR4Locus <- function(loc.pos, chromosome){
    if(loc.pos["status"] == "gain"){
      segInLoc <- segList[as.character(segList[, "chrom"]) ==
                          as.character(chromosome) & as.numeric(
                          segList[, "loc.start"]) >= as.numeric(
                          loc.pos["loc.start"]) & as.numeric(
                          segList[, "loc.end"]) <= as.numeric(
                          loc.pos["loc.end"]) & as.numeric(
                          segList[, "seg.mean"]) >= 0,
                          , drop = FALSE]
    }else{
      segInLoc <- segList[as.character(segList[, "chrom"]) ==
                          as.character(chromosome) & as.numeric(
                          segList[, "loc.start"]) >= as.numeric(
                          loc.pos["loc.start"]) & as.numeric(
                          segList[, "loc.end"]) <= as.numeric(
                          loc.pos["loc.end"]) & as.numeric(
                          segList[, "seg.mean"]) < 0,
                          , drop = FALSE]
    }
    if(nrow(segInLoc) > 0){
      mcrByLocus <- findMCR(segInLoc)
      mcrs <<- rbind(mcrs, cbind(chromosome, loc.pos["span"],
                                 loc.pos["status"], loc.pos["loc.start"],
                                 loc.pos["loc.end"], mcrByLocus))
    }
  }
  # Find altered regions based on the percentile values defined
  segList <- altered(cghmcr) <- getAlteredSegs(cghmcr)
  if(nrow(segList) == 0){
    cat(paste("cghMCR settings do not allow for the detection of any",
        "altered segments.\n"))
    cat("You may change the settings and try again.\n")
    return(NA)
  }
  # Join altered regions that are separated by less than 500 kb
  spans <- as.data.frame(mergeSegs(altered(cghmcr),
                                           gapAllowed(cghmcr)))
  mcrs <- NULL
  temp <- split.data.frame(altered(cghmcr),
                           factor(altered(cghmcr)[, "chrom"]))
  #locus <- lapply(temp, getLocus, overlap = overlap)
  for(chrom in as.vector(spans[, "chrom"])){
    #locusByChrom <- locus[[chrom]]
    junk <- apply(spans[spans[, "chrom"] == chrom, , drop = FALSE], 1,
                  getMCR4Locus, chromosome = chrom) 
  }
  colnames(mcrs) <- c("chromosome", "locus", "status", "loc.start",
                      "loc.end", "mcr", "mcr.start", "mcr.end", "samples",
                      "counts")
  rownames(mcrs) <- mcrs[, 1]
  mcrs <- mcrs[as.numeric(mcrs[, "counts"])/(ncol(DNAData(cghmcr)) - 3) >=
               recurrence(cghmcr)/100, -ncol(mcrs), drop = FALSE]
  
  return(mcrs)
}

# Takes the segment data for a givan chromosome and returns a data frame
# defining the loci identified. Segments are considered to belong to the
# same locus unless they overlap by the number of bases specified by
# variable 'overlap'.
getLocus <- function(segData, overlap = 500){
  locus <- NULL
  listBySign <- list(gain = segData[segData[, "seg.mean"] >= 0,,
                     drop = FALSE], loss = segData[segData[, "seg.mean"] < 0,
                                      , drop = FALSE])
  for(name in names(listBySign)){
    if(nrow(listBySign[[name]]) > 0){
      locus <- rbind(locus, cbind(findLocus(listBySign[[name]],
                                 overlap = overlap), status = name))
    }
  }
  return(locus)
}
                  

# Filters out segments based on given thresholds
getAlteredSegs <- function(cghmcr){
  segList <- split.data.frame(DNASeg(cghmcr), factor(DNASeg(cghmcr)[, "ID"]))
  toReturn <- NULL
  for(sample in names(segList)){
    thresholds <- quantile(DNAData(cghmcr)[, sample],
                      prob = c(alteredLow(cghmcr), alteredHigh(cghmcr)))
    temp <- segList[[sample]]
    toReturn <- rbind(toReturn, temp[temp[, "seg.mean"] < thresholds[1] |
                              temp[, "seg.mean"] > thresholds[2], ,
                              drop = FALSE])
  }
  return(toReturn)
}


mergeSegs <- function(segs, gapAllowed){
  if(nrow(segs) < 2){
    span <- c(chrom = segs[, "chrom"], status = ifelse(segs[, "seg.mean"]
              >= 0, "gain", "loss"), span = "span.1",
              loc.start = segs[,"loc.start"], loc.end = segs["loc.end"])
  }else{
    span <- NULL
    temp <- splitSegments(segs)
    for(chrom in names(temp)){
      for(sign in names(temp[[chrom]])){
        segBySign <- temp[[chrom]][[sign]]
        if(nrow(segBySign) != 0){
          span  <- rbind(span, findSpan(segBySign, gapAllowed))
        }
      }
    }
  }
  return(span)
}

# finds the segment spans on a given chromosome
findSpan <- function(segData, gapAllowed = 50000){
  span <- NULL

  chrom <- as.vector(segData[1, "chrom"])
  status <-ifelse(segData[1, "seg.mean"] >= 0, "gain", "loss")
  begin <- min(segData[, "loc.start"])
  end <- max(segData[segData[, "loc.start"] == begin, "loc.end"])
  if(end >= max(segData[, "loc.end"])){
    span <- rbind(span, c(chrom, status, "span.1", begin, end))
  }else{
    spanCounter <- 1
    while(TRUE){
      # Subsetting segments that overlap with the locus currently defined
      tempData <- segData[segData[, "loc.start"] != begin, , drop = FALSE]
      tempData <- tempData[tempData[, "loc.end"] > end &
                        tempData[, "loc.start"] <=  end + gapAllowed,
                       , drop = FALSE]
      if(nrow(tempData) < 1){
        # Segments break, define one locus
        span <- rbind(span, c(chrom, status, paste("span.", spanCounter,
                                             sep = ""), begin, end))
        spanCounter <- spanCounter + 1
        # Drop the ones have been processed
        segData <- segData[segData[, "loc.start"] > end + gapAllowed, ,
                           drop = FALSE]
        if(nrow(segData) < 1){
          break
        }else{
          # Define the begining of another locus
          begin <- min(segData[, "loc.start"])
          end <- max(segData[segData[, "loc.start"] == begin, "loc.end"])
        }
      }else{
        end <- max(tempData[, "loc.end"])
        # End of the data
        if(end >= max(segData[, "loc.end"])){
          span <- rbind(span, c(chrom, status, paste("span.", spanCounter,
                                              sep = ""), begin, end))
          break
        }
      }
    }
  }
  colnames(span) <- c("chrom", "status", "span", "loc.start", "loc.end")
  return(span)
}


splitSegments <- function(segs){
  temp <- split.data.frame(segs, factor(segs[, "chrom"]))
  return(lapply(temp, splitSegByMean))
}

splitSegByMean <- function(segs){
  posMeans <- segs[as.numeric(segs[, "seg.mean"]) >= 0, , drop = FALSE]
  negMeans <- segs[as.numeric(segs[, "seg.mean"]) < 0, , drop = FALSE]
  return(list(gain = posMeans[order(posMeans[, "loc.start"]),, drop = FALSE],
              loss = negMeans[order(negMeans[, "loc.start"]), , drop = FALSE]))
}
 

findLocus <- function(segData, overlap = 500){
  locus <- NULL
  
  begin <- min(segData[, "loc.start"])
  end <- max(segData[segData[, "loc.start"] == begin, "loc.end"])
  if(end >= max(segData[, "loc.end"])){
    locus <- rbind(locus, c("locus.1", begin, end))
  }else{
    locusCounter <- 1
    while(TRUE){
      # Subsetting segments that overlap with the locus currently defined
      tempData <- segData[segData[, "loc.start"] != begin, , drop = FALSE]
      tempData <- tempData[tempData[, "loc.end"] > end &
                        tempData[, "loc.start"] - overlap <=  end,
                       , drop = FALSE]
      if(nrow(tempData) < 1){
        # Segments break, define one locus
        locus <- rbind(locus, c(paste("locus.", locusCounter, sep = ""),
                                    begin, end))
        locusCounter <- locusCounter + 1
        # Drop the ones have been processed
        segData <- segData[segData[, "loc.start"] - overlap > end, ,
                           drop = FALSE]
        if(nrow(segData) < 1){
          break
        }else{
          # Define the begining of another locus
          begin <- min(segData[, "loc.start"])
          end <- max(segData[segData[, "loc.start"] == begin, "loc.end"])
        }
      }else{
        end <- max(tempData[, "loc.end"])
        # End of the data
        if(end >= max(segData[, "loc.end"])){
          locus <- rbind(locus, c(paste("locus.", locusCounter, sep = ""),
                                  begin, end))
          break
        }
      }
    }
  }
  colnames(locus) <- c("locus", "loc.start", "loc.end")
  return(locus)
}


# Takes the segment data for a given locus and returns the MCRs located
# by mapping the data to the potential MCRs identified by getPotentialMCR
findMCR <- function(segData){
  mcrs <- NULL
  if(nrow(segData) == 1){
    mcrs <- matrix(c("mcr.1", segData[, "loc.start"], segData[, "loc.end"],
              segData[, "ID"], 1), ncol = 5)
  }else{
    mcrCounter <- 1
    pMCRs <- getPotentialMCR(segData)
    for(index in 1:nrow(pMCRs)){
      segs <- segData[as.numeric(segData[, "loc.start"])
                      <= as.numeric(pMCRs[index, "mcr.start"]) &
                      as.numeric(segData[, "loc.end"])
                      >= as.numeric(pMCRs[index, "mcr.end"]), "ID"]
      if(length(segs) > 0){
        mcrs <- rbind(mcrs, c(paste("mcr", mcrCounter, sep = "."),
                              pMCRs[index, "mcr.start"],
                              pMCRs[index, "mcr.end"],
                              paste(unique(segs), sep = "", collapse = ","),
                              length(unique(segs))))
        mcrCounter <- mcrCounter + 1
      }
    }
  }
  colnames(mcrs) <- c("mcr.num", "mcr.start", "mcr.end", "samples", "counts")

  return(mcrs)
}


# Locate all the potential MCRs by aligning the starting and ending points
# of all the segments from all the samples. The region between any two
# adjacent points will be a potential MCR
getPotentialMCR <- function(segData){
  temp <- sort(unique(as.numeric(c(segData[, "loc.start"],
                                   segData[, "loc.end"]))))
  mcrs <- cbind(temp[1:(length(temp) - 1)], temp[2:length(temp)])
  colnames(mcrs) <- c("mcr.start", "mcr.end")

  return(mcrs)
}


mergeMCRProbes <- function(mcr, rawData){
  temp <- apply(mcr, 1, FUN = function(x){
    probes <- rawData[rawData[, 2] == x[1] & as.numeric(rawData[, 3])
                        >= as.numeric(x[7]) & as.numeric(rawData[, 3])
                        <= as.numeric(x[8]), 1]
    if(length(probes) == 0){
      return(NA)
    }else{
      return(paste(probes, sep = "", collapse = ","))
    }                   
  })

  mcr <- cbind(mcr, probes = temp)
  return(mcr)
}


plot.DNAcopy <- function(x, ..., save = FALSE, layout){
  args <- list(...)
  sampleNames <- unique(x[["output"]][, "ID"])
  locations <- cghMCR:::alignGenes(x[["data"]][, c("chrom", "maploc")])
  adjustments <- cghMCR:::getAdjustments(x[["data"]][, c("chrom", "maploc")])
  segdata <- cghMCR:::adjustSegments(x[["output"]], adjustments)
  if(save){
    tempPng <- file.path(tempdir(), "segments.png")
    png(filename = tempPng)
  }
  if(missing(layout)){
    par(mfrow = c(length(sampleNames), 1))
  }else{
    par(mfrow = layout)
  }
  for(sampleName in sampleNames){
    plot(0, 0, cex = 0, main = sampleName, xlab = "Chromsome",
         ylab = " Log2 ratio", ylim = c(-5, 5), axes = FALSE,
         xlim = c(0, max(locations) + 10))
    axis(2)
    box()
    cghMCR:::highlightChrom(adjustments)
    points(locations, x[["data"]][, sampleName], cex = 0.4, pch = 16)
    lines(c(min(locations), max(locations)), rep(0, 2), lwd = 2)
    lines(c(min(locations), max(locations)), rep(1, 2))
    lines(c(min(locations), max(locations)), rep(-1, 2))
    cghMCR:::drawSegs(segdata[segdata[, "ID"] == sampleName, ])
    cghMCR:::markChrom(adjustments)
  }
  if(save){
    dev.off()
    return(tempPng)
  }else{
    return(invisible())
  }
}

plot.MCR <- function(x, ..., DNAData, threshold = 1, save = FALSE,
                     expand = c(2000, 2000), layout, range = 1:5){
  
  x <- x[unlist(lapply(strsplit(x[, "samples"], ","), length)) >= threshold, ]

  if(save){
    tempPng <- paste(tempfile("mcrs"), "png", sep = ".")
    png(filename = tempPng)
  }
  if(!missing(layout)){
    par(mfrow = layout)
  }
  for(index in range){
    showMCR(x[index, "mcr.start"], x[index, "mcr.end"],
            DNAData[DNAData[, "chrom"] == x[index, "chromosome"] &
                     as.numeric(DNAData[, "maploc"]) >=
                     (as.numeric(x[index, "mcr.start"]) - expand[1]) &
                     as.numeric(DNAData[, "maploc"]) <=
                     (as.numeric(x[index, "mcr.end"]) + expand[2]),
                     c(3, which(colnames(DNAData) %in%
                                unlist(strsplit(x[index, "samples"], ",")))),
                     drop = FALSE], "median")
  }
  if(save){
    dev.off()
    return(tempPng)
  }else{
    return(invisible())
  }
}
  

showMCR <- function(start, end, ratioMat, what = c("mean", "median")){
  what = match.arg(what)
  temp <- cbind(ratioMat[, 1], unlist(apply(ratioMat,
                                     1, FUN = function(ratios){
                                     ratios <- ratios[2:length(ratios)]
                                     ratios<- ratios[!is.na(ratios)]
                                     if(what == "mean"){
                                       return(mean(as.numeric(ratios)))
                                     }else{
                                       return(median(as.numeric(ratios)))
                                     }
                                   })))
  left <- as.numeric(unlist(ratioMat[as.numeric(ratioMat[, 1]) <
                                     as.numeric(start),
                                     2:ncol(ratioMat), drop = FALSE]))
  inSeg <- as.numeric(unlist(ratioMat[as.numeric(ratioMat[, 1]) >=
                                      as.numeric(start) &
                                      as.numeric(ratioMat[, 1]) <=
                                      as.numeric(end),
                                      2:ncol(ratioMat), drop = FALSE]))
  right <- as.numeric(unlist(ratioMat[as.numeric(ratioMat[, 1]) >
                                      as.numeric(end),
                                      2:ncol(ratioMat), drop = FALSE]))
  plot(temp[, 1], temp[, 2], xlab = "Chromosomal location",
       ylab = "Log2 ratio")
  if(length(left) > 0){
    lines(c(min(ratioMat[as.numeric(ratioMat[, 1]) < as.numeric(start), 1]),
            max(ratioMat[as.numeric(ratioMat[, 1]) < as.numeric(start), 1])),
          rep(ifelse(what == "mean", mean(left), median(left)), 2),
          col = "red", lwd = 2)
  }
  lines(c(start, end), rep(ifelse(what == "mean", mean(inSeg), median(inSeg)),
                           2), col = "red", lwd = 2)
  if(length(right) > 0){
    lines(c(min(ratioMat[as.numeric(ratioMat[, 1]) > as.numeric(end), 1]),
            max(ratioMat[as.numeric(ratioMat[, 1]) > as.numeric(end), 1])),
          rep(ifelse(what == "mean", mean(left), median(right)), 2),
          col = "red", lwd = 2)
  }
  

  return(invisible())
}

getLineData <- function(start, end, ratioMat, what){
  left <- as.numeric(unlist(ratioMat[as.numeric(ratioMat[, 1]) < start,
                                     2:ncol(ratioMat), drop = FALSE]))
  inSeg <- as.numeric(unlist(ratioMat[as.numeric(ratioMat[, 1]) >= start &
                                      as.numeric(ratioMat[, 1]) <= end,
                                      2:ncol(ratioMat), drop = FALSE]))
  right <- as.numeric(unlist(ratioMat[as.numeric(ratioMat[, 1]) > end,
                                      2:ncol(ratioMat), drop = FALSE]))
  if(length(left) > 0){
    means <- c(means, ifelse(what == "mean", mean(left), median(left)))
  }
  means <- c(means, ifelse(what == "mean", mean(inSeg), median(inSeg)))
  if(length(right) > 0){
    means <- c(means, ifelse(what == "mean", mean(right), median(right)))
  }
  
  return(means)
}

markChrom <- function(adjustments){
  chromLocs <- NULL
  chromNames <- NULL
  for(i in 1:length(adjustments) - 1){
    if(i %% 2 == 1){
      chromLocs <- c(chromLocs, mean(c(adjustments[i], adjustments[i + 1])))
      chromNames <- c(chromNames, names(adjustments)[i])
    }
  }
  text(chromLocs, rep(-5.125, length(chromLocs)), chromNames, cex = 0.5)
  return(invisible())
}


highlightChrom <- function(adjustments){
  for(index in 1:length(adjustments)){
    if(index %% 2 == 1){
      polygon(c(adjustments[index], adjustments[index + 1],
                adjustments[index + 1], adjustments[index]),
              c(-5, -5, 5, 5), col = "gray", border = "white")
    }
  }
  return(invisible())
}


drawSegs <- function(segdata){
  drawSegLine <- function(segLocs){
    lines(c(segLocs["loc.start"], segLocs["loc.end"]),
          rep(segLocs["seg.mean"], 2), col = "red", lwd = 2)
  }
  junck <- apply(segdata, 1, FUN = drawSegLine)
  return(invisible())
}

# This function gets the values that can be used to adjust chromosomal
# locations to align genes on different chromosomes so that they  appear
# in sequence from chromosome one to Y along a single strand
getAdjustments <- function(positions){
  temp <- split.data.frame(positions, factor(positions[, 1]))
  temp <- unlist(lapply(temp, FUN = function(x) max(x[, 2])))
  chroms <- sort(as.numeric(names(temp)[names(temp) %in% 1:24]))
  if(any(names(temp) %in% "X")){
    chroms <- c(chroms, "X")
  }
  if(any(names(temp) %in% "Y")){
    chroms <- c(chroms, "Y")
  }
  adjustments <- 0
  for(index in 1:length(chroms) - 1){
    adjustments <- c(adjustments, (adjustments[length(adjustments)]+
                                   temp[as.character(chroms[index])]))
  }
  names(adjustments) <- chroms
  return(adjustments)
}

alignGenes <- function(positions){
  adjustments <- getAdjustments(positions)
  for(chrom in names(adjustments)){
    positions[positions[, 1] == chrom, 2] <-
      positions[positions[, 1] == chrom, 2] + adjustments[chrom]
  }
  return(positions[, 2])
}

adjustSegments <- function(segdata, adjustments){
  for(chrom in names(adjustments)){
    segdata[segdata[, "chrom"] == chrom, "loc.start"] <-
      segdata[segdata[, "chrom"] == chrom, "loc.start"] + adjustments[chrom]
    segdata[segdata[, "chrom"] == chrom, "loc.end"] <-
      segdata[segdata[, "chrom"] == chrom, "loc.end"] + adjustments[chrom]
  }
  return(segdata)
}


  
getSegData <- function(arrayRaw){
  filter <- getProbeFilter(arrayRaw)
  log2Ratio <- maM(arrayRaw)[filter, ]
  
  chrom <- gsub("chr([0-9XY]+):.*", "\\1", maInfo(maGnames(
                               arrayRaw))["SystematicName"][filter, 1])
  pos <- gsub("chr.*-([0-9]+)", "\\1", maInfo(maGnames(
                     arrayRaw))["SystematicName"][filter, 1])
  probes <- maInfo(maGnames(arrayRaw))["ProbeName"][filter, 1]
  CNA.object <- CNA(matrix(as.numeric(log2Ratio),
                    ncol = ncol(log2Ratio), byrow = FALSE), chrom,
                    as.numeric(pos), data.type = "logratio",
                    sampleid = colnames(log2Ratio))
  segdata <- segment(smooth.CNA(CNA.object))
  segdata[["data"]] <- cbind(probe = probes, as.data.frame(segdata[["data"]]))

  return(segdata)
}

getProbeFilter <- function(arrayRaw){
  chrom <- gsub("chr([0-9XY]+):.*", "\\1", 
                          maInfo(maGnames(arrayRaw))["SystematicName"][, 1])
  
  return(intersect(which(chrom %in% c(1:22, "X", "Y")), which(maInfo(
                   maGnames(arrayRaw))["ControlType"] >= 0)))
}


cghMCR <- function(segments, gapAllowed = 500, alteredLow = 0.03,
                   alteredHigh = 0.97, spanLimit = 20000000,
                   recurrence = 75){
  return(new("cghMCR", DNASeg = segments[["output"]], gapAllowed = gapAllowed,
             DNAData = segments[["data"]], alteredLow = alteredLow,
             alteredHigh = alteredHigh, spanLimit = spanLimit, 
             recurrence = recurrence))
}


showSegment <- function(dataMat, chromosome, start, end, samples = NA,
                        what = c("mean", "median")){
  getMean <- function(x){
    x <- x[-(1:3)]
    if(!is.na(samples)){
      x <- x[samples]
    }
    return(ifelse(what == "mean", mean(as.numeric(x), na.rm = TRUE),
                  median(as.numeric(x), na.rm = TRUE)))
  }
  what = match.arg(what)
  temp <- dataMat[dataMat[, "chrom"] == chromosome &
                  as.numeric(dataMat[, "maploc"]) >= as.numeric(start) &
                  as.numeric(dataMat[, "maploc"]) <= as.numeric(end),
                  , drop = FALSE]
  y <- apply(temp, 1, FUN = getMean)
  plot(temp[, "maploc"], y, xlab = "Chromosomal location",
       ylab = ifelse(what == "mean", "Mean across samples",
         "Median across samples"))
  
  return(invisible())
}
