# Keep identifiers and study-design fields exact while allowing harmless
# floating-point rounding in explicitly named, derived summary columns.
ccle_same_ordered_frame <- function(actual, expected,
                                    derived_columns = character(),
                                    tolerance = 1e-12) {
  if (!is.data.frame(actual) || !is.data.frame(expected) ||
      !identical(names(actual), names(expected)) ||
      nrow(actual) != nrow(expected) ||
      !all(derived_columns %in% names(expected))) return(FALSE)
  exact_columns <- setdiff(names(expected), derived_columns)
  if (!isTRUE(all.equal(
    actual[exact_columns], expected[exact_columns],
    tolerance = 0, check.attributes = FALSE
  ))) return(FALSE)
  all(vapply(derived_columns, function(column) {
    left <- actual[[column]]
    right <- expected[[column]]
    is.numeric(left) && is.numeric(right) &&
      all(is.finite(left)) && all(is.finite(right)) &&
      all(abs(left - right) <= tolerance)
  }, logical(1L)))
}

# Record the inputs, analysis sources and runtime of a newly computed artifact.
ccle_artifact_provenance <- function(input_files, code_files) {
  checksum <- function(paths) {
    if (any(!file.exists(paths))) {
      stop("Cannot record provenance for missing file(s): ",
           paste(paths[!file.exists(paths)], collapse = ", "), call. = FALSE)
    }
    values <- unname(tools::md5sum(paths))
    names(values) <- basename(paths)
    values
  }
  list(
    format_version = 1L,
    hash_algorithm = "md5",
    input_md5 = checksum(input_files),
    code_md5 = checksum(code_files),
    session_info = utils::sessionInfo()
  )
}
