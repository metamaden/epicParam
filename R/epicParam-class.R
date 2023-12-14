#!/usr/bin/env R

# Author: Sean Maden

#' epicParam-class
#'
#' Applies the EPIC::EPIC() deconvolution algorithm.
#' 
#' @details Main constructor for class \linkS4class{epicParam}.
#' @rdname epicParam-class
#' @seealso 
#' \linkS4class{referencebasedParam-class}
#' 
#' @references 
#' Racle, Julien, and David Gfeller. 2020. “EPIC: A Tool to Estimate the 
#' Proportions of Different Cell Types from Bulk Gene Expression Data.” In 
#' Bioinformatics for Cancer Immunotherapy: Methods and Protocols, edited by 
#' Sebastian Boegel, 233–48. Methods in Molecular Biology. New York, NY: 
#' Springer US. https://doi.org/10.1007/978-1-0716-0327-7_17.
#' 
#' @examples
#' exampleDataList <- getDeconvolutionExampleData()
#' newParam <- epicParam(exampleDataList[["bulkExpression"]],exampleDataList[["referenceExpression"]],exampleDataList[["cellScaleFactors"]])
#' deconvolution(newParam)
#' 
setClass("epicParam", contains="referencebasedParam", 
         slots=c(z.var = "matrix"))

#' Make new object of class epicParam
#'
#' Main constructor for class \linkS4class{epicParam}.
#'
#' @param bulkExpression Bulk mixed signals matrix of samples, which can be 
#' matched to single-cell samples.
#' @param referenceExpression Signature matrix of cell type-specific signals. 
#' If not provided, can be computed from a
#' provided ExpressionSet containing single-cell data.
#' @param cellScaleFactors Cell size factor transformations of length equal to 
#' the K cell types to deconvolve.
#' @param z.var Variance table for the reference matrix.
#' @param returnInfo Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @references 
#' Racle, Julien, and David Gfeller. 2020. “EPIC: A Tool to Estimate the 
#' Proportions of Different Cell Types from Bulk Gene Expression Data.” In 
#' Bioinformatics for Cancer Immunotherapy: Methods and Protocols, edited by 
#' Sebastian Boegel, 233–48. Methods in Molecular Biology. New York, NY: 
#' Springer US. https://doi.org/10.1007/978-1-0716-0327-7_17.
#' 
#' @export
epicParam <- function(bulkExpression, referenceExpression, 
                      cellScaleFactors=NULL, z.var=NULL, returnInfo=FALSE) {
  if(is(z.var, "NULL")){z.var <- matrix()}
  new("epicParam", bulkExpression = bulkExpression, 
      referenceExpression = referenceExpression, 
      cellScaleFactors = cellScaleFactors, z.var = z.var, 
      returnInfo = returnInfo)
}



#' Deconvolution method for class \linkS4class{epicParam}
#' 
#' Main deconvolution method for the \linkS4class{epicParam} to run the 
#' \code{EPIC::EPIC()} implementation of the EPIC algorithm.
#' 
#' @param object An object of class \linkS4class{epicParam}.
#' 
#' @examples
#' exampleDataList <- getDeconvolutionExampleData()
#' newParam <- epicParam(exampleDataList[["bulkExpression"]],exampleDataList[["referenceExpression"]],exampleDataList[["cellScaleFactors"]])
#' deconvolution(newParam)
#' 
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#'
#' @references 
#' Racle, Julien, and David Gfeller. 2020. “EPIC: A Tool to Estimate the 
#' Proportions of Different Cell Types from Bulk Gene Expression Data.” In 
#' Bioinformatics for Cancer Immunotherapy: Methods and Protocols, edited by 
#' Sebastian Boegel, 233–48. Methods in Molecular Biology. New York, NY: 
#' Springer US. https://doi.org/10.1007/978-1-0716-0327-7_17.
#'
#' @export
setMethod("deconvolution", signature(object = "epicParam"), function(object){
  require(EPIC)
  lparam <- callNextMethod()
  # instantiate and format objects
  bulkExpression <- lparam[["bulkExpression"]]
  referenceExpression <- lparam[["referenceExpression"]]
  cellScaleFactors <- lparam[["cellScaleFactors"]]
  cellScaleFactors <- as.numeric(cellScaleFactors)
  if(!"otherCells" %in% names(cellScaleFactors)){
    message("Setting size/mRNA for missing label 'otherCells' to 0...")
    cellScaleFactors["otherCells"] <- 0
  }
  referenceExpression <- as.matrix(referenceExpression)
  bulkExpression <- as.data.frame(bulkExpression)
  # get reference object
  z.var <- object[["z.var"]];z.var <- as.matrix(z.var)
  zVarCond <- ncol(z.var)==1 & nrow(z.var)==1
  if(zVarCond){
    reference <- list(refProfiles = referenceExpression,
                      sigGenes = rownames(referenceExpression))
  } else{
    reference <- list(refProfiles = referenceExpression,
                      sigGenes = rownames(referenceExpression),
                      refProfiles.var = z.var)
  }
  # get predictions
  result <- lapply(seq(ncol(bulkExpression)), function(bulkSampleIndex){
    EPIC::EPIC(bulk = bulkExpression[,bulkSampleIndex,drop=FALSE], 
               reference = reference,
               mRNA_cell = cellScaleFactors)
  })
  names(result) <- colnames(bulkExpression)
  # get standardized return list
  predictions <- lapply(result, function(iter){iter$mRNAProportions})
  returnList <- parseDeconvolutionPredictionsResults(
    predictions, c(colnames(referenceExpression), "otherCells"), 
    colnames(bulkExpression))
  if(object[["returnInfo"]]){
    returnList <- list(
      predictions=predictions,
      result.info = result,
      metadata = parametersList[["metadata"]])}
  return(returnList)
})
