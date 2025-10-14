# --- Internal helpers (same file, not exported) ------------------------------
# ATTENZIONE: STAI USANDO VARIANTI TREEROWZERO: check_matrix_criteria (influsice su agas_fitness) , agas_mutation

# --- Internal state (no globals leaked) --------------------------------------
.pkg_state <- new.env(parent = emptyenv())
.pkg_state$best_individuals_all <- list()
.pkg_state$best_individual      <- NULL
.pkg_state$best_fitness         <- -Inf


# --- Definition of constant necessary for create_sem_model_string_from_matrix 

.variables <- c("eta1", "eta2", "eta3", "eta4", "eta5", "eta6") 
.n_variables <- length(.variables)
.measurement_model <- list(
  eta1 = c("y1", "y2", "y3"),
  eta2 = c("y4", "y5", "y6"),
  eta3 = c("y7", "y8", "y9"), 
  eta4 = c("y10", "y11", "y12"),
  eta5 = c("y13", "y14", "y15"),
  eta6 = c("y16", "y17", "y18")
)
.type_of_variable <- c(eta1 = "composite", eta2 = "composite", 
                      eta3 = "composite", eta4 = "composite",
                      eta5 = "composite", eta6 = "composite")
.structural_coefficients <- list()

# --- Mutazione --------------------------------------
.agas_mutation <- function(object, parent) {
  mutate <- parent <- as.vector(object@population[parent,])
  mutate_matrix <- matrix(mutate, nrow = .n_variables, byrow = TRUE)
  
  diag(mutate_matrix) <- 0  
  
  # Ensure the first three rows are all zeros
  mutate_matrix[1, ] <- 0
  mutate_matrix[2, ] <- 0
  mutate_matrix[3, ] <- 0 
  
  mutate_vector <- as.vector(t(mutate_matrix))
  
  # Create indices of lower triangular part (excluding diagonal) starting from row 4
  indices <- which(!diag(.n_variables), arr.ind = TRUE)
  indices <- indices[indices[, 1] > 3, ]  # Exclude first three rows
  
  # Convert row and column indices to vector indices
  if (length(indices) > 0) {
    subdiag_indices <- (indices[, 1] - 1) * .n_variables + indices[, 2]
    # Select a random index from the sub-diagonal indices 
    if (length(subdiag_indices) > 0 && runif(1) <= 0.2) {
      available_indices <- subdiag_indices  # Keep track of available indices to flip
      while (length(available_indices) > 1) {
        j <- sample(available_indices, size = 1)
        mutate_vector[j] <- abs(mutate_vector[j] - 1)
        
        # Check if the mutation creates a cycle
        adj_matrix <- matrix(mutate_vector, nrow = .n_variables, byrow = TRUE)
        g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "directed", diag = FALSE)
        has_cycle <- has_cycle_dfs(g, adj_matrix)
        
        if (has_cycle) {
          mutate_vector[j] <- abs(mutate_vector[j] - 1)
          available_indices <- setdiff(available_indices, j)
        } else {
          break  
        }
      }
    }
  }
  return(mutate_vector)
}

# --- DFS cycle routine -------------------------------------------------------
#' @keywords internal
has_cycle_dfs <- function(graph, adj_matrix) {
  visited  <- rep(FALSE, igraph::vcount(graph))
  recStack <- rep(FALSE, igraph::vcount(graph))
  
  for (v in seq_len(igraph::vcount(graph))) {
    if (!visited[v]) {
      if (dfs_util(graph, v, visited, recStack, adj_matrix)) {
        return(TRUE)
      }
    }
  }
  FALSE
}

#' @keywords internal
dfs_util <- function(graph, v, visited, recStack, adj_matrix) {
  visited[v]  <- TRUE
  recStack[v] <- TRUE
  
  neighbors <- which(adj_matrix[v, ] == 1)
  for (u in neighbors) {
    if (!visited[u] && dfs_util(graph, u, visited, recStack, adj_matrix)) {
      return(TRUE)
    } else if (recStack[u]) {
      return(TRUE)
    }
  }
  
  recStack[v] <- FALSE
  FALSE
}

# --- Matrix utils ------------------------------

#' @keywords internal
.matrix_to_string <- function(m) {
  row_strs <- apply(m, 1, paste, collapse = ",")
  paste(row_strs, collapse = ";")
}

#' @keywords internal
.check_matrix_criteria <- function(adj_matrix) {
  if (any(rowSums(adj_matrix[1:3, , drop = FALSE]) != 0)) return(FALSE) # first row all zeros
  if (any(diag(adj_matrix) != 0)) return(FALSE)                       # diag all zeros
  TRUE
}

#' @keywords internal
.repair_individual_unused <- function(adj_matrix) {
  n <- nrow(adj_matrix)
  row_sums <- rowSums(adj_matrix); col_sums <- colSums(adj_matrix)
  empty_indices <- which(row_sums == 0 & col_sums == 0)
  if (length(empty_indices) > 0) {
    for (k in empty_indices) {
      if (runif(1) < 0.5) {
        valid_rows <- if (k == n) 2:(n-1) else (k+1):n 
        if (k <= n) {
          valid_rows <- valid_rows[valid_rows != k]
        }
        if (length(valid_rows) > 0) {
          if (length(valid_rows) == 1) {
            adj_matrix[valid_rows, k] <- 1
          } else {
            condition_met <- TRUE
            while (length(valid_rows) > 0 && condition_met) {
              chosen_row <- sample(valid_rows, 1)
              adj_matrix_mod <- adj_matrix
              adj_matrix_mod[chosen_row, k] <- 1
              valid_rows <- valid_rows[-which(valid_rows == chosen_row)]
              g <- igraph::graph_from_adjacency_matrix(adj_matrix_mod, mode = "directed", diag = FALSE)
              condition_met <- has_cycle_dfs(g, adj_matrix_mod)
              if (!condition_met) {
                adj_matrix <- adj_matrix_mod
                break  
              }
            }
          }
        }
      } else {
        valid_cols <- if (k == n) 1:(n-1) else (k+1):n 
        if (k <= n) {
          valid_cols <- valid_cols[valid_cols != k]
        }
        if (length(valid_cols) > 0) {
          if (length(valid_cols) == 1) {
            adj_matrix[k, valid_cols] <- 1
          } else {
            condition_met <- TRUE
            while (length(valid_cols) > 0 && condition_met) {
              chosen_col <- sample(valid_cols, 1)
              adj_matrix_mod <- adj_matrix
              adj_matrix_mod[chosen_col, k] <- 1
              valid_cols <- valid_cols[-which(valid_cols == chosen_col)]
              g <- igraph::graph_from_adjacency_matrix(adj_matrix_mod, mode = "directed", diag = FALSE)
              condition_met <- has_cycle_dfs(g, adj_matrix_mod)
              if (!condition_met) {
                adj_matrix <- adj_matrix_mod
                break  
              }
            }
          }
        }
      }
    }
  }
  return(adj_matrix)
}

#' @keywords internal
.check_matrix <- function(mat) {
  n <- nrow(mat)
  for (i in seq_len(n)) if (all(mat[i, ] == 0) && all(mat[, i] == 0)) return(FALSE)
  TRUE
}

# --- Create model string from matrix -----------------------------------------

#' @keywords internal
.create_sem_model_string_from_matrix <- function(adj_matrix) {
  model_string <- "# Composite model\n" # Initialize the model string
  
  # Include the measurement model specified in the input for composite types
  for (var in names(.measurement_model)) {
    if (.type_of_variable[var] == "composite") {
      items <- paste(.measurement_model[[var]], collapse = " + ")
      model_string <- paste(model_string, sprintf("  %s <~ %s\n", var, items), sep = "")
    }
  }
  
  # Adding reflective measurement models
  model_string <- paste(model_string, "\n# Reflective measurement model\n")
  for (var in names(.measurement_model)) {
    if (.type_of_variable[var] == "reflective") {
      items <- paste(.measurement_model[[var]], collapse = " + ")
      model_string <- paste(model_string, sprintf("  %s =~ %s\n", var, items), sep = "")
    }
  }
  
  # Structural model section
  model_string <- paste(model_string, "\n# Structural model\n")
  for (i in seq_along(.variables)) {
    dependent <- .variables[i]
    predictors <- .variables[adj_matrix[i, ] == 1]
    if (length(predictors) > 0) {
      relationship_str <- paste(predictors, collapse = " + ")
      model_string <- paste(model_string, sprintf("  %s ~ %s\n", dependent, relationship_str), sep = "")
    }
  }
  
  model_string <- gsub("\\n\\s+$", "", model_string) # Cleanup the string
  model_string <- trimws(model_string)  
  
  return(model_string)
}

# --- Fitness ---------------------------------------------------

#' @keywords internal
.bic_fitness <- function(adj_matrix, dataset_generated) {
  model_string <- .create_sem_model_string_from_matrix(adj_matrix)
  out <- csem(.data = dataset_generated, .model = model_string)
  model_criteria <- calculateModelSelectionCriteria(
    out, .by_equation = FALSE, .only_structural = FALSE
  )
  model_criteria$BIC
}

#' @keywords internal
.agas_fitness <- function(matrix_vector, dataset_generated) {
  n_variables <- length(.variables)
  adj_matrix  <- matrix(matrix_vector, nrow = n_variables, byrow = TRUE)
  
  if (!.check_matrix(adj_matrix)) {
    adj_matrix <- .repair_individual_unused(adj_matrix)
  }
  
  if (!.check_matrix_criteria(adj_matrix)) return(-100000)
  
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "directed", diag = FALSE)
  if (has_cycle_dfs(g, adj_matrix)) return(-100000)
  
  model_string <- .create_sem_model_string_from_matrix(adj_matrix)
  out <- csem(.data = dataset_generated, .model = model_string)
  ver <- verify(out)
  if (!sum(ver) == 0) return(-100000)
  
  sem_fitness <- -.bic_fitness(adj_matrix, dataset_generated)
  if (is.na(sem_fitness)) return(-100000)
  
  if (sem_fitness > .pkg_state$best_fitness) {
    .pkg_state$best_individual      <- adj_matrix
    .pkg_state$best_fitness         <- sem_fitness
    .pkg_state$best_individuals_all <- append(.pkg_state$best_individuals_all, list(adj_matrix))
  }
  
  sem_fitness
}

# --- Public function ---------------------------------------------------------

#' Do a model search
#'
#' \lifecycle{stable}
#'
#' Perform a model search \insertCite{Hair2016;textual}{cSEM}
#'
#' @usage doModelSearch(.object = NULL)
#'
#' @return A numeric value (demo returns BIC and fitness while you experiment).
#' @inheritParams csem_arguments
#' @seealso [cSEMResults]
#' @references \insertAllCited{}
#' @export
doModelSearch <- function(.object = NULL, 
                          .popsize = 100, 
                          .maxiter=50,
                          seeds = c(2, 4, 5)) {
  if (is.null(.object)) stop("`.object` must be a cSEM results object.")
  
  # compute criteria once
  model_criteria <- calculateModelSelectionCriteria(
    .object,
    .by_equation = FALSE,
    .only_structural = FALSE
  )
  
  BIC = model_criteria$BIC
  print(BIC)
  
  agas_dataset <- .object$Information$Data
  
  seeds <- as.integer(seeds)
  mat_list <- vector("list", length(seeds))
  names(mat_list) <- as.character(seeds)
  
  for (i in seq_along(seeds)) {
    ga_control <- GA::ga(
        type = "binary",
        nBits = .n_variables * .n_variables,
        popSize = .popsize,
        maxiter = .maxiter,
        pmutation = 1.0,
        pcrossover = 0.8,
        fitness = function(x) .agas_fitness(x, agas_dataset),
        elitism = TRUE,
        parallel = FALSE,
        seed = i,
        mutation = .agas_mutation
      )
    mat_list[[i]] <- .pkg_state$best_individual
  }
  
  arr <- array(unlist(mat_list),
               dim = c(.n_variables, .n_variables, length(mat_list)))
  mean_mat <- apply(arr, c(1, 2), mean, na.rm = TRUE)
  print(mean_mat)
  out <- .matrix_to_string(mean_mat)
  return(out)
  
}
