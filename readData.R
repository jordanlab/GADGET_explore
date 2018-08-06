readGrs = function(fileName, sampleInfo){
  # clean up the comment lines
  grs = my.read.table(fileName, comment.char = "##")
  colnames(grs) = c("IND", "GRS")
  grs = merge(grs, sampleInfo, by="IND")
}

my.read.table = function(fileName, comment.char){
  con = file(fileName)
  open(con)
  clean.lines = sub(paste0(comment.char, ".*"), "", readLines(con))
  close(con)
  clean.lines = sub("SAMPLE.*", "", clean.lines)
  read.table(text = paste(clean.lines, collapse = "\n"))
}