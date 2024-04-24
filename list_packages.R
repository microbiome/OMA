# Function to read and parse the DESCRIPTION file
extract_packages_from_file <- function(file_path) {
  # Read the entire file
  lines <- readLines(file_path)
  
  # Helper function to extract packages from the 'Suggests' section
  extract_suggests_packages <- function(lines) {
    section_name <- "Suggests"
    section_index <- grep(paste0("^", section_name, ":"), lines)
    if (length(section_index) > 0) { # Check if the section exists
      start_index <- section_index + 1
      # Find the index of the next section or end of file if no next section
      end_index <- ifelse(any(grep("^[A-Za-z]+:", lines[-(1:section_index)])),
                          min(grep("^[A-Za-z]+:", lines[-(1:section_index)])) + section_index - 1,
                          length(lines))
      section_lines <- lines[start_index:end_index]
      
      # Concatenate, split by commas, and trim whitespace
      section_string <- paste(section_lines, collapse = " ")
      packages <- unlist(strsplit(section_string, ","))
      packages <- trimws(packages)  # Remove leading and trailing whitespaces
      # Filter out any empty entries
      return(packages[packages != ""])
    } else {
      return(character(0))  # Return empty if no section found
    }
  }
  
  # Extract packages from the Suggests section
  suggests_packages <- extract_suggests_packages(lines)
  
  return(suggests_packages)
}

# Path to the DESCRIPTION file
desc_file_path <- "DESCRIPTION"

# Extract packages and save to CSV
packages <- extract_packages_from_file(desc_file_path)
packages_df <- data.frame(Package = packages,stringsAsFactors = FALSE)
write.table(packages, "oma_packages.csv", row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)

# Output message indicating completion
cat("Packages from the 'Suggests' section have been extracted and saved to 'oma_packages.csv'.\n")
