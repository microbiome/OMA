extract_packages <- function(description_file_path = "~/OMA/DESCRIPTION", output_file_path = "~/OMA/oma_packages/oma_packages.csv") {
  # Check if the DESCRIPTION file exists
  if (!file.exists(description_file_path)) {
    stop(paste("DESCRIPTION file not found at path:", description_file_path))
  }
  
  # Read the DESCRIPTION file
  lines <- tryCatch(readLines(description_file_path), error = function(e) {
    message("Error reading DESCRIPTION file: ", e$message)
    return(NULL)
  })
  
  if (is.null(lines)) return()  # Exit if file could not be read
  
  # Function to extract packages from a given section
  extract_packages_from_section <- function(lines, section_name) {
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
      message("Section '", section_name, "' not found in DESCRIPTION file.")
      return(character(0))
    }
  }
  
  # Extract packages from the Imports and Suggests sections
  imports_packages <- extract_packages_from_section(lines, "Imports")
  suggests_packages <- extract_packages_from_section(lines, "Suggests")
  
  # Combine the packages into a single vector and sort alphabetically
  combined_packages <- sort(c(imports_packages, suggests_packages))
  
  # Write the combined packages to a CSV file without headers
  tryCatch({
    write.table(combined_packages, output_file_path, row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
    message("Combined packages written to ", output_file_path)
  }, error = function(e) {
    message("Error writing combined packages to CSV file: ", e$message)
  })
}

# Call the function with the default path to the DESCRIPTION file and output file path
extract_packages()
