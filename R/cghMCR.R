# A class that contains the segmentation and MCR derived from the
# segmentation data
# 
# Copyright 2006, Jianhua Zhang, all rights reserved
#

setClass("cghMCR", representation(DNASeg = "data.frame",
                                  DNAData = "data.frame",
                                  altered = "data.frame",
                                  gapAllowed = "numeric",
                                  spanLimit = "numeric",
                                  alteredLow = "numeric",
                                  alteredHigh = "numeric",
                                  recurrence = "numeric",
                                  thresholdType = "character"))

setGeneric("DNASeg",
           function(object) standardGeneric("DNASeg"))
setMethod("DNASeg", "cghMCR",
          function(object) object@DNASeg)

setGeneric("thresholdType",
           function(object) standardGeneric("thresholdType"))
setMethod("thresholdType", "cghMCR",
          function(object) object@thresholdType)

setGeneric("DNAData",
             function(object) standardGeneric("DNAData"))
setMethod("DNAData", "cghMCR",
          function(object) object@DNAData)

setGeneric("altered",
           function(object) standardGeneric("altered"))
setMethod("altered", "cghMCR",
          function(object) object@altered)

setGeneric("gapAllowed",
           function(object) standardGeneric("gapAllowed"))
setMethod("gapAllowed", "cghMCR",
          function(object) object@gapAllowed)

setGeneric("alteredLow",
           function(object) standardGeneric("alteredLow"))
setMethod("alteredLow", "cghMCR",
          function(object) object@alteredLow)

setGeneric("alteredHigh",
           function(object) standardGeneric("alteredHigh"))
setMethod("alteredHigh", "cghMCR",
          function(object) object@alteredHigh)

setGeneric("spanLimit",
           function(object) standardGeneric("spanLimit"))
setMethod("spanLimit", "cghMCR",
          function(object) object@spanLimit)

setGeneric("recurrence",
           function(object) standardGeneric("recurrence"))
setMethod("recurrence", "cghMCR",
          function(object) object@recurrence)

setGeneric("spans",
           function(object) standardGeneric("spans"))
setMethod("spans", "cghMCR",
          function(object) object@spans)
  
setMethod("show", "cghMCR",
          function(object) {
            cat("Object of cghMCR\n")
            cat("DNASeg:\n")
            cat(paste("\tRow = ", nrow(DNASeg(object)), "\n", sep = ""))
            cat(paste("\tColumn = ", ncol(DNASeg(object)), "\n", sep = ""))
            if(nrow(DNASeg(object)) <= 5){
              rows <- nrow(DNASeg(object))
            }else{
              rows <- 5
            }
            cat(DNASeg(object)[1:rows, ])
            if(nrow(DNASeg(object) > 5)){
              cat("\n ..........\n")
            }
            cat("Parameter settings:\n")
            for(fun in c("gapAllowed", "spanLimit", "alteredLow",
                         "alteredHigh", "recurrence")){ 
              cat(paste("\t", fun, " = ", do.call(fun, args = list(object)),
                        "\n", sep = ""))
            }
          })

setGeneric("MCR",
           function(object) standardGeneric("MCR"))
setMethod("MCR", "cghMCR",
          function(object) {
            temp <- getMCR(object)
            class(temp) <- "MCR"
            return(temp)
          })

####### Replace methods
setGeneric("DNASeg<-", function(object, value)
           standardGeneric("DNASeg<-"))
setReplaceMethod("DNASeg", "cghMCR", function(object, value){
  object@DNASeg <- value; object})

setGeneric("DNAData<-", function(object, value)
           standardGeneric("DNAData<-"))
setReplaceMethod("DNAData", "cghMCR", function(object, value){
  object@DNAData <- value; object})

setGeneric("thresholdType<-", function(object, value)
           standardGeneric("thresholdType<-"))
setReplaceMethod("thresholdType", "cghMCR", function(object, value){
  object@thresholdType <- value; object})

setGeneric("altered<-", function(object, value)
           standardGeneric("altered<-"))
setReplaceMethod("altered", "cghMCR", function(object, value){
  object@altered <- value; object})

setGeneric("gapAllowed<-", function(object, value)
           standardGeneric("gapAllowed<-"))
setReplaceMethod("gapAllowed", "cghMCR", function(object, value){
  object@gapAllowed <- value; object})

setGeneric("alteredLow<-", function(object, value)
           standardGeneric("alteredLow<-"))
setReplaceMethod("alteredLow", "cghMCR", function(object, value){
  object@alteredLow <- value; object})

setGeneric("alteredHigh<-", function(object, value)
           standardGeneric("alteredHigh<-"))
setReplaceMethod("alteredHigh", "cghMCR", function(object, value){
  object@alteredHigh <- value; object})

setGeneric("recurrence<-", function(object, value)
           standardGeneric("recurrence<-"))
setReplaceMethod("recurrence", "cghMCR", function(object, value){
  object@recurrence <- value; object})

setGeneric("spanLimit<-", function(object, value)
           standardGeneric("spanLimit<-"))
setReplaceMethod("spanLimit", "cghMCR", function(object, value){
  object@spanLimit <- value; object})

setGeneric("spans<-", function(object, value)
           standardGeneric("spans<-"))
setReplaceMethod("spans", "cghMCR", function(object, value){
  object@spans <- value; object})


# A method for segmenting arrayCGH data added for marrayRaw
#setGeneric("getSegments",
#           function(object) standardGeneric("getSegments"))
#setMethod("getSegments", "marrayRaw",
#          function(object) getSegData(object))

#setMethod("getSegments", "marrayNorm",
#          function(object) getSegData(object))


# Class and methods for Segment Gain Or Loss
setClass("SGOL", representation(gol = "data.frame",
                                threshold = "vector",
                                method = "function"))

setGeneric("gol",
           function(object) standardGeneric("gol"))
setMethod("gol", "SGOL",
          function(object) object@gol)
setGeneric("gol<-", function(object, value)
           standardGeneric("gol<-"))
setReplaceMethod("gol", "SGOL", function(object, value){
  object@gol <- value; object})

setGeneric("threshold",
           function(object) standardGeneric("threshold"))
setMethod("threshold", "SGOL",
          function(object) object@threshold)
setGeneric("method",
           function(object) standardGeneric("method"))
setMethod("method", "SGOL",
          function(object) object@method)
          
setMethod("[", "SGOL", function(x, i, j, ..., drop = FALSE) {
    if (missing(drop)) drop <- FALSE
    if (missing(i) && missing(j)) {
        if (length(list(...))!=0)
            stop("specify genes or SGOL scores to subset")
          return(x)
    }
    if (!missing(j) && !missing(i)){
        gol(x) <- gol(x)[i, j, ..., drop = drop]
    }else{
        if (!missing(i)){
            gol(x) <- gol(x)[i,, ..., drop=drop]
        }else{
            gol(x) <- gol(x)[, j, ..., drop = drop]
        }
    }    
    return(x)
})

setMethod("colnames", "SGOL",
          function(x, do.NULL = TRUE, prefix = "col"){
              colnames(col(x), do.NULL = do.NULL, prefix = prefix)
})

setMethod("rownames", "SGOL",
          function(x, do.NULL = TRUE, prefix = "row"){
              rownames(col(x), do.NULL = do.NULL, prefix = prefix)
})
         
setGeneric("plot")
setMethod("plot", "SGOL",
          function(x, y, ...){
            if(!missing(y)) {
                if(is.logical(y)){
                    plotSGOL(gol(x), XY = y)
                }else{
                    plotSGOL(gol(x))
                }
            }else{
                plotSGOL(gol(x))
            }
})

setGeneric("GEOI",
           function(object) standardGeneric("GEOI"))
setMethod("GEOI", "SGOL",
          function(object) getGEOI(gol(object)))

