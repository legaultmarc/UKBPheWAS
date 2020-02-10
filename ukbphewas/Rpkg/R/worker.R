library(rzmq)
library(rjson)


read_and_cast_feather <- function(filename) {
  df <- as.data.frame(arrow::read_feather(filename))

  for (i in 1:ncol(df)) {
    if (names(df)[i] != "sample_id") {
      df[, i] <- as.numeric(df[, i])
    }
  }

  df
}


get_xs <- function(configuration) {
  get_reader <- function(filename) {
    if (endsWith(filename, ".feather")) {
      return(read_and_cast_feather)
    } else if (endsWith(filename, ".csv") || endsWith(filename, ".csv.gz")) {
      return(read.csv)
    } else {
      stop("Invalid covariable file.")
    }
  }

  # rjson parses a list with a single string as a character vector.
  # We have to explicitly handle it.
  if (length(configuration$covars_filenames) == 1) {
    filename <- configuration$covars_filenames
    data <- get_reader(filename)(filename)
    data$sample_id <- as.character(data$sample_id)
    return(data)
  }

  # This is the case where multiple data files need to be joined together.
  out <- NULL
  for (filename in configuration$covars_filenames) {
    cat(paste0("R: Reading covars ('", filename, "').\n"))

    cur <- get_reader(filename)(filename)
    cur$sample_id <- as.character(cur$sample_id)

    if (is.null(out)) {
      out <- cur
    } else {
      out <- merge(out, cur, by.x = "sample_id", by.y = "sample_id",
                   all.x = T, all.y = T)
    }
  }

  out$sample_id <- as.character(out$sample_id)

  out
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
