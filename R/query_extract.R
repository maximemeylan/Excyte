query_extract <- function(processed_fcs_obj,channels=c("all","with_desc")[1]){
  all_channels <- processed_fcs_obj$all_channels
  processed_fcs_df <- processed_fcs_obj$processed_fcs
  #if channels are not specified
  if (all(channels == "with_desc")) {
    channels_to_select <- all_channels[!is.na(all_channels[, "desc"]), 1]
    processed_fcs_df <- processed_fcs_df[, c(channels_to_select, "sample_id")]
    channels <- colnames(processed_fcs_df)[colnames(processed_fcs_df) != "sample_id"]
  }
  if (all(channels != "with_desc") & all(channels != "all")) {
    processed_fcs_df <- processed_fcs_df[, c(channels, "sample_id")]
    channels <- colnames(processed_fcs_df)[colnames(processed_fcs_df) != "sample_id"]
  } else{
    channels <- all_channels[,1]
  }
  processed_fcs_df <- processed_fcs_df[,c(channels,"sample_id")]
  attr(processed_fcs_df,"all_channels") <- all_channels
  return(processed_fcs_df)
}
