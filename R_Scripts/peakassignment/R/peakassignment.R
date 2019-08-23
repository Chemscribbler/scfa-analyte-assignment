#' Assign Analytes
#'
#' Works through a join (on ColID), so scaling may be an issue with large datasets
#'
#' @param rttable the table of retention time values
#' @param peaktable the table of peak values
#'
#' @return Returns the analyte (as a factor) or an NA value
#' @export
assign_analytes <- function(peaktable, rttable){
  if (!"Exp" %in% colnames(peaktable)) {
    peaktable <- sample_label(peaktable)
  }
  output <- left_join(rttable, peaktable, by = c("ColID")) %>%
    filter(RTmin < RT, RTmax > RT) %>%
    select(ColID, Analyte, SID, Exp, Plate, Well,
           Inst, InjTime, RT, PeakWidth, Area,
           SampleType, CurveID, PlateID, ControlNum) %>%
    filter(!is.na(Analyte))
  return(output)
}

#' Close misses checker
#'
#' Meant to check for peaks that are just outside of the window
#'
#' @param df raw peak data, no manipulation
#' @param rttable this will be constructed from the raw data if not submitted
#' (should be equaivalent to calling the function yourself)
#' @param comparisontable If you want to call out a specific dataset, otherwise
#' it will build it from the raw submitted data equivalently to the standards method
#'
#' @return returns the anti-join table, which (mathematically) should be just peaks
#' that are not labelled in the dataset
#' @export
#'
check_nearby_peaks <- function(df, rttable=NA, comparisontable=NA){
  if(is.na(rttable)){
    rttable <- df %>% generate_rt_table()
  }
  if (is.na(comparisontable)) {
    comparisontable <- df %>% assign_analytes(rttable)
  }
  for (ColID in unique(rttable$ColID)) {
    for (AnalyteID in unique(rttable$Analyte)) {
      rttable <- adjust_rt(rttable, ColID, AnalyteID,adjustMax = FALSE, -0.02)
      rttable <- adjust_rt(rttable, ColID, AnalyteID,adjustMax = TRUE, 0.02)
    }
  }
  temp <- assign_analytes(df, rttable)
  return(anti_join(comparisontable,temp))
}
