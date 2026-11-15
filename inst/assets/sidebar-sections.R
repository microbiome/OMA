# Post-render step: list each chapter's sections (e.g. 3.1, 3.2) in the left
# sidebar, under the chapter link, on every page of the book.
#
# Quarto's book sidebar can only link whole pages, so after rendering we read
# the numbered level-2 headings from each HTML page and insert them as links
# into the sidebar of every page. The script can be run several times: blocks
# inserted by a previous run are removed first.

out_dir <- Sys.getenv("QUARTO_PROJECT_OUTPUT_DIR", "docs")
block_class <- "sidebar-chapter-sections"

read_html <- function(f) {
    x <- readChar(f, file.size(f), useBytes = TRUE)
    Encoding(x) <- "UTF-8"
    x
}

# Numbered level-2 headings of one page: number, anchor id and title
get_sections <- function(html) {
    pattern <- paste0(
        '(?s)<h2 data-number="([^"]+)"[^>]*data-anchor-id="([^"]+)"[^>]*>\\s*',
        '<span class="header-section-number">[^<]*</span>\\s*(.*?)\\s*</h2>'
    )
    m <- regmatches(html, gregexpr(pattern, html, perl = TRUE))[[1]]
    if (length(m) == 0) return(NULL)
    parts <- regmatches(m, regexec(pattern, m, perl = TRUE))
    data.frame(
        number = vapply(parts, `[`, "", 2),
        id = vapply(parts, `[`, "", 3),
        # Keep inline code, drop other tags such as footnote links
        title = gsub("<(?!/?code\\b)[^>]*>", "",
            vapply(parts, `[`, "", 4), perl = TRUE),
        stringsAsFactors = FALSE
    )
}

sections_html <- function(sec, href) {
    items <- sprintf(
        '<li><a href="%s#%s"><span class="section-number">%s</span> %s</a></li>',
        href, sec$id, sec$number, sec$title
    )
    paste0('<ul class="', block_class, '">', paste(items, collapse = ""),
        "</ul>")
}

add_sections <- function() {
    files <- list.files(out_dir, pattern = "\\.html$", recursive = TRUE)
    pages <- files[startsWith(files, "pages/")]

    # Collect the sections of each chapter page
    sections <- list()
    for (p in pages) {
        sec <- get_sections(read_html(file.path(out_dir, p)))
        if (!is.null(sec)) sections[[p]] <- sec
    }

    # Chapter links in the sidebar, e.g. ../pages/containers.html
    link_pattern <- paste0(
        '(?s)<a href="([^"]*?)(pages/[^"#]+\\.html)" ',
        'class="sidebar-item-text sidebar-link[^"]*">.*?</a>\\s*</div>'
    )
    old_block <- paste0('(?s)<ul class="', block_class, '">.*?</ul>')

    for (f in files) {
        path <- file.path(out_dir, f)
        html <- read_html(path)
        if (!grepl('id="quarto-sidebar"', html, fixed = TRUE)) next
        html <- gsub(old_block, "", html, perl = TRUE)

        hits <- gregexpr(link_pattern, html, perl = TRUE)
        links <- regmatches(html, hits)[[1]]
        if (length(links) == 0) next
        new <- vapply(links, function(link) {
            parts <- regmatches(link, regexec(link_pattern, link, perl = TRUE))[[1]]
            sec <- sections[[parts[3]]]
            if (is.null(sec)) return(link)
            paste0(link, sections_html(sec, paste0(parts[2], parts[3])))
        }, "", USE.NAMES = FALSE)
        regmatches(html, hits) <- list(new)

        con <- file(path, open = "wb")
        writeChar(enc2utf8(html), con, eos = NULL, useBytes = TRUE)
        close(con)
    }
}

# Do not fail the book build if this optional step breaks
tryCatch(add_sections(), error = function(e) {
    warning("Could not add chapter sections to the sidebar: ",
        conditionMessage(e), call. = FALSE)
})
