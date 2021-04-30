library(rzmq)
library(rjson)


# Get the covariables from the databundle.
#
# This function looks at a data source called 'worker_data'.
# It also applies subsetting if required.
get_xs <- function(configuration) {
  cat(paste(
    "R: Reading worker data from databundle '", configuration$databundle_path,
    "'\n"
  ))
  data <- databundle.deserialize.parquet(
    configuration$databundle_path, "worker_data"
  )

  if (!is.null(configuration$subset)) {
    cat(paste0(
        "R: Subsetting ", length(configuration$subset), " individuals.\n"
    ))
    data <- subset(data, sample_id %in% configuration$subset)
  }

  return(data)
}


remove_sex_from_rhs <- function(rhs, sex_column) {
  rhs_li <- char_to_terms(rhs)

  rhs_li <- rhs_li[rhs_li != sex_column]
  paste(rhs_li, collapse = " + ")
}


char_to_terms <- function(rhs) {
  labels(terms(as.formula(paste0("~", rhs))))
}


# Remove columns that have no variance.
#
# This is typically useful if, for example, sex is included in the right hand
# side but individuals of one sex get excluded by the data generator.
#
# The freeze.col.idx argument is used so that some columns (e.g. the first one
# corresponding to the outcome) do not get affected by the filter.
drop_columns_with_no_variance <- function(mat,
                                          freeze.col.idx = 1,
                                          verbose = T) {
  n_unique <- apply(mat, 2, function(x) length(unique(x)))

  keep_cols <- n_unique > 1

  # Apply the freeze
  keep_cols[freeze.col.idx] <- TRUE
  dropped_cols <- colnames(mat)[!keep_cols]

  if (verbose && length(dropped_cols) > 0) {
    cat(paste0("R: Dropping columns with no variance: '", dropped_cols, "'\n"))
  }

  list(
    dropped_cols = dropped_cols,
    mat = mat[, keep_cols]
  )
}


Worker <- function(worker_id, dealer_addr, monitor_addr, callback) {
  # Connect to the dealer.
  context <- init.context()
  socket <- init.socket(context, "ZMQ_DEALER")
  connect.socket(socket, dealer_addr)
  cat(paste0("R: Connected to dealer socket on '", dealer_addr, "'\n"))

  # Connect to the monitor and send ready signal.
  monitor_sock <- init.socket(context, "ZMQ_REQ")
  connect.socket(monitor_sock, monitor_addr)
  send.socket(monitor_sock,
              charToRaw(paste0("ready ", worker_id)),
              serialize=F)

  res <- receive.socket(monitor_sock, unserialize=F)
  res <- rawToChar(res)
  stopifnot(res == "OK")

  # Start the main worker loop.
  while(T) {

    res <- receive.socket(socket, unserialize=F, dont.wait=T)

    if (is.null(res)) {

      # Ask the monitor if the main process is done pushing data.
      send.socket(monitor_sock, charToRaw("is_done?"), serialize=F)
      is_done <- receive.socket(monitor_sock, unserialize=F)

      if (length(is_done) == 0) {
        # We need to wait for data (not done).
        Sys.sleep(1)
        next
      } else {
        cat("R: no more data to process, exiting.\n")
        break
      }

    }

    else {
      # A multipart message containing metadata and then data should have
      # been sent, we now collect the parts.
      parts <- list(res)

      while(get.rcvmore(socket)) {
        parts <- append(parts, list(receive.socket(socket, unserialize=F)))
      }

      stopifnot(length(parts) == 2)

      # These messages should be metadata and then data.
      metadata <- fromJSON(rawToChar(parts[[1]]))

      # Time callback.
      t0 <- proc.time()
      callback(metadata, parts[[2]])
      elapsed <- proc.time() - t0
      print(elapsed)
    }

  }

  warnings()
  cat(paste0("R: worker ", worker_id, " done.\n"))
}
