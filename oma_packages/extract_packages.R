packages <- function() {
  lines <- readLines("../DESCRIPTION")
  
  extract_packages <- function(lines, section_name) {
    section_index <- grep(paste0("^", section_name, ":"), lines)
    if (length(section_index) > 0) {
      start_index <- section_index + 1
      end_index <- ifelse(any(grep("^[A-Za-z]+:", lines[-(1:section_index)])),
                          min(grep("^[A-Za-z]+:", lines[-(1:section_index)])) + section_index - 1,
                          length(lines))
      section_lines <- lines[start_index:end_index]
      section_string <- paste(section_lines, collapse = " ")
      packages <- unlist(strsplit(section_string, ","))
      packages <- trimws(packages)
      return(packages[packages != ""])
    } else {
      return(character(0))
    }
  }
  
  suggests_packages <- extract_packages(lines, "Suggests")
  write.table(suggests_packages, "oma_packages.csv", row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
}

packages()
