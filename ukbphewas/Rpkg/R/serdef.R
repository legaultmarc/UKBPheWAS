library(rjson)


dtypeToRType <- list(
    float64 = list(type="double", size=8),
    float32 = list(type="double", size=4),
    int64 = list(type="integer", size=8),
    int32 = list(type="integer", size=4)
)


# This is the deserialization function from the zeromq queue.
# TODO: We could use a better representation from arrow here instead.
deserialize <- function(data) {
    # First uint is metadata length.
    metadata_length <- readBin(data, what="integer", n=1)
    metadata <- fromJSON(rawToChar(data[5:(4+metadata_length)]))

    data <- data[(5+metadata_length):length(data)]
    
    li <- list()
    for (meta in metadata) {
        cur <- data[1:(meta$n_bytes)]
        data <- data[(meta$n_bytes+1):length(data)]

        if (meta$dtype == "string") {
            # Parse string
            bytes_per_str <- meta$n_bytes / meta$n

            v <- rep("", times=meta$n)
            idx <- 1
            for (i in seq(1, meta$n, bytes_per_str)) {
              v[idx] <- rawToChar(cur[i:(i+bytes_per_str-1)])
              idx <- idx + 1
            }

        } else {
            # Parse numerical data
            # Get the right type metadata.
            t <- dtypeToRType[[meta$dtype]]
            v <- readBin(cur, t$type, n=meta$n, size=t$size)
        }

        li[[meta$name]] <- v
    }
    
    as.data.frame(li, stringsAsFactors=FALSE)
}


# See https://github.com/legaultmarc/databundle
# for main project description and up to date bindings.
databundle.deserialize.parquet <- function(databundle_filename, data_source_name) {

  if (!endsWith(data_source_name, ".parquet")) {
    data_source_name <- paste0(data_source_name, ".parquet")
  }

  # Untar relevant member.
  tmpdir_name <- tempdir()
  untar(databundle_filename, files=data_source_name, exdir=tmpdir_name)
  parquet_filename <- file.path(tmpdir_name, data_source_name)

  data <- arrow::read_parquet(parquet_filename)

  file.remove(parquet_filename)

  return(data)

}
