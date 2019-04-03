#' Replace data in Spectrum
#' 
#' Replace the spectral data and peak annotations in a Spectrum object with a dataframe.
#' Sets all slots which are present as columns in the given dataframe.
#' 
#' Note for the behavior when mismatching columns are present:
#' 
#' If \code{clean=TRUE}, columns missing in \code{df} but present in \code{peakAnnotations}
#' are deleted from \code{o}. If \code{addNew=TRUE}, columns present in \code{df} and missing
#' in \code{peakAnnotations} are automatically created with the class of the column.
#' Otherwise, \code{strict=TRUE} raises an error; \code{strict=FALSE} ignores new and 
#' deleted columns. Strict mode is mandatory if the number of peaks changes.
#' 
#' @name setData
#' @aliases setData,Spectrum,data.frame-method
#' 
#' @param o The \code{Spectrum} object to modify
#' @param df The data frame with new data
#' @param clean \code{TRUE} if annotation columns which aren't present as 
#'  columns in the data frame should be cleared.
#' @param strict If \code{TRUE}, no annotations can be silently dropped from \code{df}
#'  or silently ignored
#' @return The modified \code{Spectrum}.
#' 
#' @author stravsmi
#' @docType methods
#' @export
setMethod("setData", c("Spectrum", "data.frame"), function(o, df, 
  clean = FALSE, addNew = FALSE, strict=TRUE)
{
  if(!strict & (nrow(df) != o@peaksCount))
    stop("If peaks are added or removed, strict mode is mandatory")
  o <- .setData.slots(o, df, clean, strict)
  o <- .setData.annotations(o, df, clean, addNew, strict)
  o
})

.setData.slots <- function(o, df, clean = FALSE, strict=TRUE)
{
  # Which columns from the data frame go to slots?
  cols <- c("mz" = "mz", "i" = "intensity" )
  o@peaksCount <- as.integer(nrow(df))
  
  replacement <- colnames(df)
  slotsMatching <- intersect(names(cols), replacement)
  slotsMissing <- setdiff(names(cols), replacement)
  
  for(col in slotsMatching)
  {
    if(!is.numeric(df[,col]))
      stop(paste0("Incorrect data type in replacement data: ", col))
    slot(o, cols[col]) <- df[,col]
  }
  
  if(length(slotsMissing) > 0)
  {
    if(strict & !clean)
      stop("mz or i is missing in replacement data. Use strict=FALSE or clean=TRUE if desired")
    if(clean)
      for(col in slotsMissing)
      {
        slot(o, col) <- rep(NA_real_, o@peaksCount)
      }
  }
  o
}

.setData.annotations <- function(o, df, clean = FALSE, addNew = FALSE, strict=TRUE)
{
  # first set  mz, i
  # then find which columns have corresponding annotations
  df$mz <- NULL
  df$i <- NULL
  
  # Find all annotation columns which have a column in the new data (df),
  # and see how to replace
  existingAnn <- colnames(o@peakAnnotations)
  replacementAnn <- colnames(df)
  
  matchingAnn <- intersect(existingAnn, replacementAnn)
  deletedAnn <- setdiff(existingAnn, replacementAnn)
  newAnn <- setdiff(replacementAnn, existingAnn)
  
  # For the matching annotations, check class then set
  for(ann in matchingAnn)
  {
    if(class(df[,ann]) != class(o@peakAnnotations[,ann]))
      stop(paste0("Annotation class mismatch for ", ann, ": ",
                  class(df[,ann]),class(o@peakAnnotations[,ann])))
      o@peakAnnotations[,ann] <- df[,ann]
  }
  
  # For the new annotations:
  # * in strict mode: fail if addNew is disallowed
  # * in loose mode: silently drop columns if addNew is disallowed
  # * if addNew: add annotation columns for new columns
  if(length(newAnn) > 0)
  {
    if(strict & !addNew)
      stop("In strict mode, no implicit dropping of peak annotation columns is permitted. Set strict=FALSE or addNew=TRUE")
    if(addNew)
      for(ann in newAnn)
      {
        # todo: improve this one for performance
        o <- addPeakAnnotation(o, newAnn, class(df[,ann]), df[,ann])
      }
  }

  # For the dropped annotations:
  # * in strict mode: fail if clean is disabled
  # * in loose mode with clean=FALSE: silently keep annotations
  # * if clean=TRUE: delete the annotation columns if they are missing in df
  if(length(deletedAnn) > 0)
  {
    if(strict & !clean)
      stop("In strict mode, no implicit preservation of existing peakAnnotation columns is permitted. Set strict=FALSE or clean=TRUE")
    if(clean)
      for(ann in deletedAnn)
        o@peakAnnotations[,ann] <- NULL
  }
  o
}


#' @export
.selectPeaks <- setMethod("selectPeaks", c("Spectrum"), function(o, filter, ..., enclos=parent.frame(2))
{
  if(missing(filter))
    return(o)

  df <- as.data.frame(o)
  f <- substitute(filter)
  df <- df[eval(f, df, enclos) & !is.na(eval(f, df, enclos)),,drop=FALSE]
  o <- setData(o, df)

  o
})
#' @describeIn addProperty Add a new column to the RmbSpectrum2 properties
#'
#' @export 
.addPeakAnnotation <- setMethod("addPeakAnnotation", 
          c("Spectrum", "character", "character", "ANY"), 
          function(o, name, type, value=NA, FUN=NULL)
{
  if(!is.null(FUN))
  {
    value <- FUN(as.data.frame(o))
  }
  if(ncol(o@peakAnnotations) == 0)
    o@peakAnnotations <- data.frame(row.names = seq_len(o@peaksCount))
  if(length(value) == 1)
    o@peakAnnotations[,name] <- as(rep(value, o@peaksCount), type)
  else if(length(value) == o@peaksCount)
    o@peakAnnotations[,name] <- value
  else
    stop("Incorrect length for annotation")
  o
})

#setGeneric("setData",	function(s, df, ...) standardGeneric("setData"))


#' @export
setMethod("peakAnnotation", c("Spectrum", "character"), function(o, name)
{
  if(peakAnnotation %in% colnames(o@peakAnnotations))
    return(o@peakAnnotation[,name])
  else
    # We can't use FALSE or NA, since it could be confused with a 1-length logical FALSE or 1-length ANY NA  
    warning("Selected an inexistent peakAnnotation column")
    return(NULL)
})


.peakAnnotationSet <- function(o, name, value, FUN = NULL, addNew = FALSE, class="")
{
  if(!is.null(FUN))
  {
    value <- FUN(as.data.frame(o))
  }
  if(class == "") class <- class(value)
  if(length(value) != o@peaksCount)
    stop("Incorrect number of entries in value")
  if(!(name %in% colnames(o@peakAnnotation)) & !addNew)
  {
    warning("Trying to set inexistent annotation. To enable autogeneration, set addNew = TRUE.")
    return(o)
  }
  else if(!(name %in% colnames(o@peakAnnotations)))
    o <- addPeakAnnotation(o, name, class)
  
  if(length(value) == 1)
    o@peakAnnotations[,name] <- as(rep(value, o@peaksCount), type)
  else if(length(value) == o@peaksCount)
    o@peakAnnotations[,name] <- value
  else
    stop("Incorrect length for annotation")
  
  return(o)
}

#' @export
setMethod("peakAnnotation<-", c("Spectrum", "character", "ANY", "logical", "character"), .peakAnnotationSet )
#' @export
setMethod("peakAnnotation<-", c("Spectrum", "character", "ANY", "missing", "character"), .peakAnnotationSet )
#' @export
setMethod("peakAnnotation<-", c("Spectrum", "character", "ANY", "logical", "missing"), .peakAnnotationSet)
#' @export
setMethod("peakAnnotation<-", c("Spectrum", "character", "ANY", "missing", "missing"), .peakAnnotationSet )

# Read the existing annotation columns and types
#' @export
setMethod("peakAnnotations", c("Spectrum"),
          function(o)
            {
            name <- colnames(o@peakAnnotations)
            class <- unlist(lapply(
              name,
              class(o@peakAnnotations[,name])
            )) 
            
            data.frame(name=name, class=class) 
            })
          
