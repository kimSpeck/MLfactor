createFolder <- function(folderPath) {
  if (!file.exists(folderPath)){
    dir.create(folderPath)
  }
}

loadRData <- function(fileName){
  ###
  # load rda file and assign it to any chosen object name instead of the original 
  #   object name
  ###
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# function to parse file name to extract condition parameters
parseFilename <- function(x) {
  # pattern captures signed integers or decimals after each tag
  pat <- paste0(
    "^.*_([^_]+)",                      
    "_miItems([+-]?\\d+)",
    "_P([+-]?\\d+)",
    "_Q([+-]?\\d+)",
    "_cAlpha([+-]?\\d*\\.?\\d+)",
    "_varT([+-]?\\d*\\.?\\d+)",
    "_corGO([+-]?\\d*\\.?\\d+)",
    "_R2([+-]?\\d*\\.?\\d+)\\.rds$"
  )
  m <- stringr::str_match(x, pat)
  if (is.na(m[1,1])) stop("Pattern not matched.")
  out <- as.data.frame(m[, -1, drop = FALSE], stringsAsFactors = FALSE)
  names(out) <- c("MIvio","miItems","P","Q","cAlpha","varT","corGO","R2")
  # convert numeric fields (leave MIvio as character)
  num_cols <- setdiff(names(out), "MIvio")
  out[num_cols] <- lapply(out[num_cols], as.numeric)
  out["MIvio"] <- substring(out["MIvio"], 6, nchar(out["MIvio"]))
  out
}


themeFunctionPaper <- function(plotObject, guides = T){
  pTmp <- plotObject + theme(
    panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
    panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.background = element_rect(color = "white", fill = "white"),
    plot.title = element_text(size = 30, face = "bold"),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    #legend.position = c(.85, .39), 
    legend.position = "bottom", 
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.key.width = unit(2, "cm"),
    legend.box = "horizontal")
  if (guides == F) {
    pTmp <- pTmp + guides(color = "none", shape = "none", linetype = "none")
  }
  return(pTmp)
}

