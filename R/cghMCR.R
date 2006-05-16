# A class that contains the segmentation and MCR derived from the
# segmentation data
# 
# Copyright 2006, Jianhua Zhang, all rights reserved
#

setClass("cghMCR", representation(DNASeg = "data.frame",
                                  DNAData = "data.frame",
                                  altered = "data.frame",
                                  spans = "data.frame",
                                  gapAllowed = "numeric",
                                  spanLimit = "numeric",
                                  alteredLow = "numeric",
                                  alteredHigh = "numeric",
                                  recurrence = "numeric"))

if(!isGeneric("DNASeg")){
  setGeneric("DNASeg",
             function(object) standardGeneric("DNASeg"))
}
setMethod("DNASeg", "cghMCR",
          function(object) object@DNASeg)

if(!isGeneric("DNAData")){
  setGeneric("DNAData",
             function(object) standardGeneric("DNAData"))
}
setMethod("DNAData", "cghMCR",
          function(object) object@DNAData)

if(!isGeneric("altered")){
  setGeneric("altered",
             function(object) standardGeneric("altered"))
}
setMethod("altered", "cghMCR",
          function(object) object@altered)

if(!isGeneric("gapAllowed")){
  setGeneric("gapAllowed",
             function(object) standardGeneric("gapAllowed"))
}
setMethod("gapAllowed", "cghMCR",
          function(object) object@gapAllowed)

if(!isGeneric("alteredLow")){
  setGeneric("alteredLow",
             function(object) standardGeneric("alteredLow"))
}
setMethod("alteredLow", "cghMCR",
          function(object) object@alteredLow)

if(!isGeneric("alteredHigh")){
  setGeneric("alteredHigh",
             function(object) standardGeneric("alteredHigh"))
}
setMethod("alteredHigh", "cghMCR",
          function(object) object@alteredHigh)

if(!isGeneric("spanLimit")){
  setGeneric("spanLimit",
             function(object) standardGeneric("spanLimit"))
}
setMethod("spanLimit", "cghMCR",
          function(object) object@spanLimit)

if(!isGeneric("recurrence")){
  setGeneric("recurrence",
             function(object) standardGeneric("recurrence"))
}
setMethod("recurrence", "cghMCR",
          function(object) object@recurrence)

if(!isGeneric("spans")){
  setGeneric("spans",
             function(object) standardGeneric("spans"))
}
setMethod("spans", "cghMCR",
          function(object) object@spans)
  
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
            for(fun in c("gapAllowed", "spanLimit", "alteredLow",
                         "alteredHigh", "recurrence")){ 
              cat(paste("\t", fun, " = ", do.call(fun, args = list(x)),
                        "\n", sep = ""))
            }
          })

if(!isGeneric("MCR")){
  setGeneric("MCR",
             function(object) standardGeneric("MCR"))
}
setMethod("MCR", "cghMCR",
          function(object) {
            temp <- getMCR(object)
            class(temp) <- "MCR"
            return(temp)
          })

####### Replace methods
if(!isGeneric("DNASeg<-")){
  setGeneric("DNASeg<-", function(object, value)
             standardGeneric("DNASeg<-"))
}
setReplaceMethod("DNASeg", "cghMCR", function(object, value){
  object@DNASeg <- value; object})

if(!isGeneric("DNAData<-")){
  setGeneric("DNAData<-", function(object, value)
             standardGeneric("DNAData<-"))
}
setReplaceMethod("DNAData", "cghMCR", function(object, value){
  object@DNAData <- value; object})

if(!isGeneric("altered<-")){
  setGeneric("altered<-", function(object, value)
             standardGeneric("altered<-"))
}
setReplaceMethod("altered", "cghMCR", function(object, value){
  object@altered <- value; object})

if(!isGeneric("gapAllowed<-")){
  setGeneric("gapAllowed<-", function(object, value)
             standardGeneric("gapAllowed<-"))
}
setReplaceMethod("gapAllowed", "cghMCR", function(object, value){
  object@gapAllowed <- value; object})

if(!isGeneric("alteredLow<-")){
  setGeneric("alteredLow<-", function(object, value)
             standardGeneric("alteredLow<-"))
}
setReplaceMethod("alteredLow", "cghMCR", function(object, value){
  object@alteredLow <- value; object})

if(!isGeneric("alteredHigh<-")){
  setGeneric("alteredHigh<-", function(object, value)
             standardGeneric("alteredHigh<-"))
}
setReplaceMethod("alteredHigh", "cghMCR", function(object, value){
  object@alteredHigh <- value; object})

if(!isGeneric("recurrence<-")){
  setGeneric("recurrence<-", function(object, value)
             standardGeneric("recurrence<-"))
}
setReplaceMethod("recurrence", "cghMCR", function(object, value){
  object@recurrence <- value; object})

if(!isGeneric("spanLimit<-")){
  setGeneric("spanLimit<-", function(object, value)
             standardGeneric("spanLimit<-"))
}
setReplaceMethod("spanLimit", "cghMCR", function(object, value){
  object@spanLimit <- value; object})

if(!isGeneric("spans<-")){
  setGeneric("spans<-", function(object, value)
             standardGeneric("spans<-"))
}
setReplaceMethod("spans", "cghMCR", function(object, value){
  object@spans <- value; object})


# A method for segmenting arrayCGH data added for marrayRaw
if(!isGeneric("getSegments")){
  setGeneric("getSegments",
             function(object) standardGeneric("getSegments"))
}
setMethod("getSegments", "marrayRaw",
          function(object) getSegData(object))

setMethod("getSegments", "marrayNorm",
          function(object) getSegData(object))

