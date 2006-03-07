# A class that contains the segmentation and MCR derived from the
# segmentation data
# 
# Copyright 2006, Jianhua Zhang, all rights reserved
#

setClass("cghMCR", representation(DNASeg = "data.frame",
                               margin = "numeric",
                               gain.threshold = "numeric",
                               loss.threshold = "numeric"))

if(!isGeneric("DNASeg")){
  setGeneric("DNASeg",
             function(object) standardGeneric("DNASeg"))
}
setMethod("DNASeg", "cghMCR",
          function(object) object@DNASeg)

if(!isGeneric("margin")){
  setGeneric("margin",
             function(object) standardGeneric("margin"))
}
setMethod("margin", "cghMCR",
          function(object) object@margin)

if(!isGeneric("gain.threshold")){
  setGeneric("gain.threshold",
             function(object) standardGeneric("gain.threshold"))
}
setMethod("gain.threshold", "cghMCR",
          function(object) object@gain.threshold)

if(!isGeneric("loss.threshold")){
  setGeneric("loss.threshold",
             function(object) standardGeneric("loss.threshold"))
}
setMethod("loss.threshold", "cghMCR",
          function(object) object@loss.threshold)


setMethod("print", "cghMCR",
          function(x, ...) {
            cat("Object of cghMCR\n")
            cat("DNASeg:\n")
            cat(paste("\tRow = ", nrow(DNASeg(x)), "\n", sep = ""))
            cat(paste("\tColumn = ", ncol(DNASeg(x)), "\n", sep = ""))
            if(nrow(DNASeg(x)) <= 5){
              rows <- nrow(DNASeg(x))
            }else{
              rows <- 5
            }
            cat(DNASeg(x)[1:rows, ])
            if(nrow(DNASeg(x) > 5)){
              cat("\n ..........\n")
            }
            cat("Parameter settings:\n")
            for(fun in c("margin", "gain.threshold", "loss.threshold")){ 
              cat(paste("\t", fun, " = ", do.call(fun, args = list(x)),
                        "\n", sep = ""))
            }
          })

if(!isGeneric("MCR")){
  setGeneric("MCR",
             function(object) standardGeneric("MCR"))
}
setMethod("MCR", "cghMCR",
          function(object) getMCR(DNASeg(object), overlap = margin(object),
                                  ampLimit = gain.threshold(object),
                                  delLimit = loss.threshold(object)))

####### Replace methods
if(!isGeneric("DNASeg<-")){
  setGeneric("DNASeg<-", function(object, value)
             standardGeneric("DNASeg<-"))
}
setReplaceMethod("DNASeg", "cghMCR", function(object, value){
  object@DNASeg <- value; object})

if(!isGeneric("margin<-")){
  setGeneric("margin<-", function(object, value)
             standardGeneric("margin<-"))
}
setReplaceMethod("margin", "cghMCR", function(object, value){
  object@margin <- value; object})

if(!isGeneric("gain.threshold<-")){
  setGeneric("gain.threshold<-", function(object, value)
             standardGeneric("gain.threshold<-"))
}
setReplaceMethod("gain.threshold", "cghMCR", function(object, value){
  object@gain.threshold <- value; object})

if(!isGeneric("loss.threshold<-")){
  setGeneric("loss.threshold<-", function(object, value)
             standardGeneric("loss.threshold<-"))
}
setReplaceMethod("loss.threshold", "cghMCR", function(object, value){
  object@loss.threshold <- value; object})


# A method for segmenting arrayCGH data added for marrayRaw
if(!isGeneric("getSegments")){
  setGeneric("getSegments",
             function(object) standardGeneric("getSegments"))
}
setMethod("getSegments", "marrayRaw",
          function(object) getSegData(object))

setMethod("getSegments", "marrayNorm",
          function(object) getSegData(object))

