#!/usr/bin/env Rscript


#This script adds genus-consistent exact-match SILVA 144 species labels to classified 16S ASVs and records a summary of the results.

suppressPackageStartupMessages({ #Loads the required packages without printing their startup messages.
  library(Biostrings) #Provides DNA FASTA reading and sequence-storage functions.
  #dada2::assignSpecies is loaded directly from its namespace only when exact matching is performed.
  library(dplyr) #Provides joins, transformations and table summaries.
  library(readr) #Provides TSV reading and writing functions.
  library(stringr) #Provides taxonomy parsing and text-matching functions.
})


#This function joins exact-match results to the classified taxonomy and adds a species label only when the classified and reference genera agree.
append_species_labels <- function(taxonomy, species_matches) {
  taxonomy %>%
    left_join(species_matches, by = "ASV_hash") %>% #Adds exact-match genus and species results to their corresponding ASVs.
    mutate(
      Classified_Genus = str_match( #Extracts the classified genus from each semicolon-delimited taxonomy string.
        Taxonomy,
        "(?:^|;\\s*)g__([^;]*)" #Captures the text following the `g__` rank prefix.
      )[, 2] %>% trimws() %>% na_if(""), #Removes surrounding whitespace and converts empty genera to missing values.
      Is_SILVA = str_detect(Taxonomy, "(?:^|;\\s*)d__"), #Identifies taxonomy strings using SILVA-style domain prefixes.
      Genus_Agrees = #Records whether the exact-match species can be accepted.
        Is_SILVA & #Requires a SILVA-formatted classification.
        !is.na(Exact_Match_Genus) & #Requires a genus from the exact sequence match.
        !is.na(Exact_Match_Species) & #Requires an unambiguous species from the exact sequence match.
        !is.na(Classified_Genus) & #Requires a genus from the classifier.
        str_to_lower(Classified_Genus) == str_to_lower(Exact_Match_Genus), #Requires the two genera to match without regard to capitalization.
      Species = if_else( #Constructs the accepted binomial species name.
        Genus_Agrees, #Adds a name only when the genus-consistency check passes.
        paste(Exact_Match_Genus, Exact_Match_Species), #Combines the exact-match genus and species epithet.
        NA_character_ #Leaves the species missing when the checks fail.
      ),
      Taxonomy = if_else( #Appends accepted species names to the taxonomy.
        !is.na(Species) & !str_detect(Taxonomy, "(?:^|;\\s*)s__"), #Requires an accepted species and avoids replacing an existing species rank.
        paste0(str_remove(Taxonomy, "\\s*;\\s*$"), "; s__", Species), #Removes a trailing separator and appends the species rank.
        Taxonomy #Otherwise preserves the original taxonomy.
      )
    )
}


#This function repeatedly attempts to copy the shared SILVA reference to node-local temporary storage for stable and faster access.
stage_reference_locally <- function(reference_path, attempts = 30L, delay = 2) {
  for (attempt in seq_len(attempts)) { #Retries staging up to the configured number of attempts.
    staged_path <- tempfile( #Creates a unique destination in the current temporary directory.
      pattern = "silva_v144_assignSpecies_", #Uses a recognizable filename prefix.
      tmpdir = tempdir(), #Places the file in node-local temporary storage when available.
      fileext = ".fa.gz" #Preserves the compressed FASTA extension.
    )
    source_size <- if (file.exists(reference_path)) { #Checks whether the shared reference is currently available.
      file.info(reference_path)$size #Records its expected size for copy validation.
    } else {
      NA_real_ #Records a missing size when the source is unavailable.
    }
    copied <- !is.na(source_size) && isTRUE(suppressWarnings( #Attempts the copy only when the source size is known.
      file.copy(reference_path, staged_path, overwrite = TRUE) #Copies the reference to the temporary location.
    ))
    staged_size <- if (file.exists(staged_path)) { #Checks whether the staged copy was created.
      file.info(staged_path)$size #Records the staged file’s size.
    } else {
      NA_real_ #Records a missing size when staging failed.
    }

    if (copied && staged_size > 0 && staged_size == source_size) { #Accepts only a nonempty copy whose size matches the source.
      message(
        "Using node-local SILVA species reference: ", staged_path, #Reports the staged reference path.
        " (", staged_size, " bytes)." #Reports the validated file size.
      )
      return(staged_path) #Returns the validated local reference path.
    }

    unlink(staged_path) #Removes an incomplete or mismatched staged copy.
    if (attempt < attempts) { #Waits and retries unless this was the final attempt.
      message(
        "Shared SILVA species reference was unavailable during staging; ",
        "retrying (", attempt, " of ", attempts, ")." #Reports the failed attempt and retry limit.
      )
      Sys.sleep(delay) #Pauses before trying the shared reference again.
    }
  }

  stop(
    "Could not make a stable node-local copy of the SILVA species reference: ",
    reference_path #Reports the reference that could not be staged.
  )
}


#This function assigns exact SILVA species matches in batches to limit peak memory use.
assign_species_in_chunks <- function(sequences, reference_path, chunk_size) {
  chunk_size <- as.integer(chunk_size) #Converts the requested batch size to an integer.
  if (is.na(chunk_size) || chunk_size < 1) { #Rejects missing, zero or negative batch sizes.
    stop("SILVA species-assignment chunk size must be a positive integer.")
  }

  sequence_count <- length(sequences) #Counts the eligible ASV sequences.
  chunk_starts <- seq.int(1L, sequence_count, by = chunk_size) #Calculates the starting index of each batch.
  chunk_results <- vector("list", length(chunk_starts)) #Preallocates storage for the batch results.

  for (chunk_number in seq_along(chunk_starts)) { #Processes each batch in order.
    first_index <- chunk_starts[[chunk_number]] #Gets the first sequence index in this batch.
    last_index <- min(first_index + chunk_size - 1L, sequence_count) #Limits the final batch to the available sequences.
    chunk <- sequences[first_index:last_index] #Selects the sequences for this batch.

    message(
      "Assigning SILVA species for ASVs ", first_index, "-", last_index,
      " of ", sequence_count,
      " (batch ", chunk_number, " of ", length(chunk_starts), ")." #Reports batch progress.
    )
    assignments <- dada2::assignSpecies( #Searches the reference for exact species-level matches.
      as.character(chunk), #Converts the DNA sequence set to character strings.
      reference_path, #Uses the node-local SILVA species reference.
      allowMultiple = FALSE, #Rejects sequences with multiple possible species assignments.
      tryRC = TRUE, #Also searches the reverse complement of each sequence.
      verbose = TRUE #Prints detailed assignment progress to the log.
    )
    chunk_results[[chunk_number]] <- tibble( #Converts this batch’s result matrix into a table.
      ASV_hash = names(chunk), #Retains the ASV identifier associated with each sequence.
      Exact_Match_Genus = assignments[, "Genus"], #Stores the exact-match reference genus.
      Exact_Match_Species = assignments[, "Species"] #Stores the exact-match species epithet.
    )

    rm(assignments, chunk) #Releases the batch data before the next iteration.
    invisible(gc(full = TRUE)) #Requests a full garbage-collection pass to reduce memory accumulation.
  }

  bind_rows(chunk_results) #Combines all batch results into one table.
}


#This function controls taxonomy loading, eligibility filtering, exact matching, genus validation and output generation.
assign_silva_species <- function(
  taxonomy_path, #Path to the exported genus-level taxonomy table.
  sequence_path, #Path to the representative-sequence FASTA.
  reference_path, #Path to the SILVA 144 assignSpecies reference.
  output_path, #Destination for the final taxonomy table.
  summary_path, #Destination for the species-assignment summary.
  chunk_size = 500L, #Maximum number of ASVs processed per reference-matching batch.
  enabled = TRUE, #Controls whether exact species matching is performed.
  silva_version = "144" #Records the selected SILVA release.
) {
  #An older classifier must never be combined with the release-144 exact-match reference.
  #When matching is disabled, the original taxonomy is copied without requiring the FASTA or reference.
  #The summary records why matching was skipped instead of reporting absent calculations as zero.
  if (!isTRUE(enabled)) { #Uses the pass-through branch when exact matching is disabled.
    if (!isTRUE(file.copy(taxonomy_path, output_path, overwrite = TRUE))) { #Copies the original taxonomy to the final destination.
      stop("Could not copy the taxonomy table while skipping species matching.")
    }
    write_tsv(
      tibble(
        status = "skipped", #Marks species assignment as intentionally skipped.
        reason = "Exact species matching is currently supported only for SILVA 144", #Explains the release restriction.
        silva_version = as.character(silva_version) #Records the configured SILVA version.
      ),
      summary_path, #Writes the skip summary to the declared output.
      na = "" #Writes missing values as empty fields.
    )
    return(invisible(NULL)) #Ends the function without loading sequences or the reference.
  }

  taxonomy <- read_tsv( #Loads the exported taxonomy table.
    taxonomy_path,
    show_col_types = FALSE, #Suppresses automatic column-type messages.
    name_repair = "minimal" #Preserves the original column names.
  )
  if (ncol(taxonomy) < 3 || !"taxonomy" %in% names(taxonomy)) { #Requires an ASV identifier, taxonomy and confidence column.
    stop("The exported taxonomy table must contain ASV ID, taxonomy, and confidence columns.")
  }
  names(taxonomy)[1] <- "ASV_hash" #Standardizes the first column as the ASV identifier.
  taxonomy <- taxonomy %>% rename(Taxonomy = taxonomy) #Renames the taxonomy column for internal processing.

  sequences <- readDNAStringSet(sequence_path) #Loads the representative sequences from FASTA.
  query_ids <- names(sequences) #Extracts the ASV identifiers from the FASTA headers.
  if (is.null(query_ids) || any(query_ids == "") || anyDuplicated(query_ids)) { #Requires every sequence to have a unique, nonempty identifier.
    stop("The representative-sequence FASTA must contain unique ASV IDs.")
  }

  classified_genera <- str_match( #Extracts each classifier-assigned genus.
    taxonomy$Taxonomy,
    "(?:^|;\\s*)g__([^;]*)" #Captures the value following the genus prefix.
  )[, 2] %>% trimws() %>% na_if("") #Cleans whitespace and converts empty genera to missing values.
  eligible_ids <- taxonomy$ASV_hash[ #Selects ASVs eligible for exact species matching.
    str_detect(taxonomy$Taxonomy, "(?:^|;\\s*)d__") & #Requires SILVA-style domain notation.
      !is.na(classified_genera) #Requires an existing classifier-assigned genus.
  ]
  eligible_sequences <- sequences[names(sequences) %in% eligible_ids] #Retains FASTA sequences belonging to eligible ASVs.

  if (length(eligible_sequences) == 0) { #Handles datasets with no eligible SILVA ASVs.
    species_matches <- tibble(
      ASV_hash = character(), #Creates an empty ASV identifier column.
      Exact_Match_Genus = character(), #Creates an empty exact-match genus column.
      Exact_Match_Species = character() #Creates an empty exact-match species column.
    )
  } else {
    local_reference_path <- stage_reference_locally(reference_path) #Copies and validates the reference in local temporary storage.
    on.exit(unlink(local_reference_path), add = TRUE) #Removes the temporary reference when the function exits.
    species_matches <- assign_species_in_chunks( #Performs batched exact matching for eligible sequences.
      eligible_sequences,
      local_reference_path,
      chunk_size
    )
  }

  augmented <- append_species_labels(taxonomy, species_matches) #Adds accepted species labels and audit columns to the taxonomy.
  summary <- augmented %>%
    summarise(
      total_16S_ASVs = n(), #Counts all ASVs in the taxonomy table.
      eligible_SILVA_ASVs = sum(Is_SILVA & !is.na(Classified_Genus)), #Counts SILVA ASVs with a classified genus.
      exact_unambiguous_species_matches = sum(!is.na(Exact_Match_Species)), #Counts unambiguous exact species matches before genus validation.
      species_labels_added = sum(Genus_Agrees, na.rm = TRUE), #Counts exact matches accepted after the genus check.
      exact_matches_withheld_by_genus_check = sum( #Counts exact matches rejected because their genus did not agree.
        !is.na(Exact_Match_Species) & !Genus_Agrees,
        na.rm = TRUE
      ),
      unmatched_or_ambiguous =
        eligible_SILVA_ASVs - exact_unambiguous_species_matches #Counts eligible ASVs without an unambiguous exact match.
    ) %>%
    mutate(status = "completed", silva_version = as.character(silva_version)) #Records successful completion and the SILVA release.

  output <- augmented %>%
    select(ASV_hash, taxonomy = Taxonomy, confidence) #Retains the ASV ID, final taxonomy and classifier confidence.
  names(output)[1] <- "#OTUID" #Restores the identifier heading expected by downstream BIOM tools.

  write_tsv(output, output_path, na = "") #Writes the final taxonomy table with blank missing values.
  write_tsv(summary, summary_path, na = "") #Writes the species-assignment audit summary.
}


#This block runs the species-assignment function when the script is executed through Snakemake.
if (exists("snakemake")) { #Prevents automatic execution when the functions are sourced outside Snakemake.
  sink(snakemake@log[[1]], append = TRUE, split = TRUE) #Sends messages to the rule log while retaining console output.
  on.exit(sink(), add = TRUE) #Restores normal output handling when execution ends.
  species_reference <- snakemake@input[["species_reference"]] #Retrieves the optional SILVA species-reference input.
  assign_silva_species(
    snakemake@input[["taxonomy"]], #Passes the exported genus-level taxonomy.
    snakemake@input[["sequences"]], #Passes the representative-sequence FASTA.
    if (length(species_reference)) species_reference[[1]] else NULL, #Passes the reference path only when it was supplied.
    snakemake@output[["taxonomy"]], #Passes the final taxonomy output path.
    snakemake@output[["summary"]], #Passes the summary output path.
    snakemake@params[["chunk_size"]], #Passes the maximum number of sequences per batch.
    enabled = snakemake@params[["enabled"]], #Passes whether exact matching is enabled.
    silva_version = snakemake@params[["silva_version"]] #Passes the configured SILVA release.
  )
}
