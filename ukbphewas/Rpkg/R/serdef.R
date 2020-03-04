library(rjson)


dtypeToRType <- list(
    float64 = list(type="double", size=8),
    float32 = list(type="double", size=4),
    int64 = list(type="integer", size=8),
    int32 = list(type="integer", size=4)
)


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
