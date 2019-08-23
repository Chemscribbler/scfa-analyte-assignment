#'Pull GC related Metadata
#'
#'`samplelabel` is a convenience function that pulls out meta-data from the
#'extracted datafile. Notably it splits up the SID column into component columns
#'of Experiment ID (EID), plate identity (PID), and Well
#'
#'Based on how the sample name is formated it will also determine if the sample
#'is one of blank/standard/calibration/sample and give the appropriate label
#'
#'
#'@param df A dataframe containing SCFA data extracted with at least the following
#'columns: SID, ColID
#'
#'
#'@return newfile a dataframe that has had metadata added to it

sample_label <- function(df){
  df %>%
    separate(col = SID,into = c("Exp","Plate","Well"),extra = "merge",fill = "right",remove = FALSE) %>%
    mutate(ColID = as.numeric(substr(ColID, 10, 11)), #Current format is ColumnID=##, so chars 10-11 corrispond to those numbers
           SampleType = if_else(Exp == "Blank", "Blank",
                                if_else(grepl("c",Plate),"Calibration",
                                        if_else(grepl("[A,C,M,T][Q,M,N,W]",Well),"Control","Sample"))),
           CurveID = if_else(grepl("Cal",SampleType),as.numeric(substr(Plate,2,3)),NaN),
           PlateID = if_else(!grepl("Calibration",SampleType),as.numeric(substr(Plate,2,3)),NaN),
           ControlNum = if_else(grepl("Con",SampleType),as.numeric(substr(Well,4,5)),NaN))
}


#' Find median retention times from standards
#'
#' @param df dataframe that contains all peaks from an SCFA GC run
#' @param analytelist a list of analytes that are in the standards, must be entered in order
#' @param excludedpoints a list of points that should be excluded
#' @param areacuttoff All points with areas less than this number will be excluded
#' @param frontrtcutoff Points with retentions times less than this number will be excluded from the front column
#' @param rearrtcutoff Points with retention times less than thsi number will be excluded from the rear column
#'
#'
#' @return a dataframe with each GC column having all analytes assigned
#'


find_median_rt <- function(df,analytelist=c("AA","PA","BA","IS"),
                                     excludedpoints=NA, areacutoff=40,
                                     frontrtcutoff=2.2, rearrtcutoff=1.8){
  scfa_levels <- c("AA","FA","PA","IsoBA","BA","IsoVA","VA","IS","IsoCA","CA","HA")
  analyte_table <- tibble(AnalyteIDNum = 1:length(analytelist),
                          Analyte= factor(analytelist, levels = scfa_levels))

  ref_table <- df %>%
    filter(!SID %in% excludedpoints,
           Area > areacutoff,
           !is.na(RT),
           if_else(ColID %% 2 == 0, RT > frontrtcutoff, RT > rearrtcutoff),
           RT < 3.8,
           SampleType == "Control") %>%
    group_by(ColID,SID) %>%
    arrange(RT, .by_group = TRUE) %>%
    mutate(AnalyteIDNum = 1:n()) %>%
    left_join(analyte_table) %>%
    ungroup() %>%
    group_by(ColID,Analyte) %>%
    summarise(
      MedRT = median(RT)
    )
  return(ref_table)
}

#' Creating Minimum and Maximum RT windows
#'
#' @param df data frame of retention times
#' @param window size of RT window, will be added and subtracted from median
#'
#' @return new reference table, will not be checkd so that a value is retuned
#'
#'
create_min_max_rt <- function(df, window=.022){
  ref_table <- df %>%
    mutate(
      RTmin = MedRT - window,
      RTmax = MedRT + window
    ) %>%
    select(ColID, Analyte, RTmin, RTmax) %>%
    filter(!is.na(Analyte))

  return(ref_table)
}

#' Adjust Retention Time
#'
#' This function can be used to adjust an individual retention time value
#'
#' @param rttable the table of retention time values
#' @param adjustColumn which ColID you want to make the adjustment on (numeric)
#' @param adjustAnalyte which analyte you want to make the adjustment on (string)
#' @param adjustMax Whether you want to adjust the MaxRT or the MinRT (True/False)
#' @param adjustment Size of the adjustment, adds in both cases, so to lower the
#' minimum RT put a negative value in.
#'
#' @return calls the check_rt_table function to verify that the table will not
#' cause other problems, and then returns either the error message or the adjusted
#' table
#'
#' @export

adjust_rt <- function(rttable, adjustColumn, adjustAnalyte,
                      adjustMax=TRUE, adjustment){
  rttable <- rttable
  if(adjustMax){
    rttable <- rttable %>%
      mutate(RTmax = if_else(ColID == adjustColumn & adjustAnalyte == Analyte,
                             RTmax+adjustment, RTmax))
  }
  else{
    rttable <- rttable %>%
      mutate(RTmin = if_else(ColID == adjustColumn & adjustAnalyte == Analyte,
                             RTmin+adjustment, RTmin))
  }
  return(check_rt_table(rttable))
}

#' Error catching method, flags potetial problems with retention time tables
#'
#'Three different conditions are checked: That the minimum is always increasing,
#'the maximumum is always inreasing, and that the maximum of an analyte is less
#'then the minmum of the following analyte
#'
#' @param df the retention table you are checking
#'
#' @export
#' @return will return the original table if there are no errors
#'

check_rt_table <- function(df){
  for (i in seq_along(df$ColID)) {
    if(!is.na(df$ColID[i+1]) && df$ColID[i] == df$ColID[i+1]){
      if (df$RTmin[i] > df$RTmin[i+1]){
        stop(paste("RTmin is greater for row ", i, "than row ", i+1))
      }
      if (df$RTmax[i] > df$RTmax[i+1]){
        stop(paste("RTmax is greater for row ", i, "than row ", i+1))
      }
      if (df$RTmax[i] > df$RTmin[i+1]) {
        stop(paste("RTmax is greater for row ", i, "than RTmin for row", i+1))
      }
    }
  }
  return(df)
}


#' Generate Retention Time Table
#'
#' This function takes raw data generated by the python script
#' and creates a retention time table by calling a few utility
#' functions. It is meant to be somewhat adjustable, and other
#' functions can be used to modify the table.
#'
#' @param df Raw data, each peak gets its own row
#' @param RTwindow the +/- from the median RT for every analyte that will be adjusted
#' @param analytelist list of analytes that are present in the controls
#' @param excludedpoints list of points (untested as of 8/21)
#' @param areacutoff will exclude all peaks below this area count
#' @param frontrtcutoff cuts off all peaks before a certain retention time on the front (even) GC column
#' @param rearrtcutoff cuts off al peaks before a certain retention time on the rear (odd) GC column
#'
#' @return returns the retention time table
#' @export
#'
#' @seealso find_median_rt
#' @seealso create_min_max_rts

generate_rt_table <- function(df, RTwindow=0.022,
                              analytelist=c("AA","PA","BA","IS"),
                              excludedpoints=NA, areacutoff=40,
                              frontrtcutoff=2.2, rearrtcutoff=1.8){
  df %>% peakassignment:::sample_label() %>%
    peakassignment:::find_median_rt(analytelist = analytelist, excludedpoints = excludedpoints,
                   areacutoff = areacutoff, frontrtcutoff = frontrtcutoff,
                   rearrtcutoff = rearrtcutoff) %>%
    peakassignment:::create_min_max_rt(window = RTwindow)
}
