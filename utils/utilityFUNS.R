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
