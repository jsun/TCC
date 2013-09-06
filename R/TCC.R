TCC <- setRefClass(
    "TCC",
    fields = list(
        count = "matrix",           # counts data of libraries.
        gene_id = "character",      # gene names
        group = "data.frame",       # groups, libraries, conditions.
        norm.factors = "numeric",   # normalization factors.
        stat = "list",              # the result of identify DE genes.
        estimatedDEG = "numeric",   # identified result by identifyDEG().
        simulation = "list",        # the aurgument inputed.
        DEGES = "list",             # detailed informations about DEGES .
        private = "list"             # intermediate data on DEGES process.
    ),

    ##  Class Methods.
    methods = list(
        initialize = function(count = NULL, group = NULL,
                              norm.factors = NULL, gene_id = NULL) {
            if (!is.null(count)) {
                ## Set up group data.
                if (is.null(group)) {
                    stop("TCC::ERROR: The group or replicates must be provided.\n")
                    #.self$group <<- data.frame(group = rep(1:length(replicates), times = replicates))
                } else {
                    if (!is.data.frame(group)) 
                        .self$group <<- data.frame(group = group)
                    else 
                        .self$group <<- group
                }
                ## Set up count data.
                if (!is.matrix(count))
                    .self$count <<- as.matrix(count)
                else
                    .self$count <<- count
                ## Set up names.
                if (is.null(rownames(.self$count))) {
                    .self$gene_id <<- paste("gene_", 
                                            c(1:nrow(count)), sep = "")
                    rownames(.self$count) <<- paste("gene_", 
                                                    c(1:nrow(count)), sep = "")
                } else {
                    .self$gene_id <<- rownames(count)
                }
                if (is.null(colnames(.self$count))) {
                    g <- as.numeric(table(group))
                    colnames(.self$count) <<- paste("G", 
                                                    rep(1:length(g), times = g),
                                                    "_rep", sequence(g), 
                                                    sep = "")
                } else {
                    colnm <- colnames(count)
                    if (sum(match(colnm, colnm)) != sum(1:length(colnm))) {
                        message("TCC::INFO: Changed the column names of count data to unique.")
                        colnames(count) <<- paste(colnm, 1:length(colnm), 
                                                  sep = "_")
                    }
                }
                rownames(.self$group) <<- colnames(.self$count)
                ## Set up normlization factors if it was given.
                if (is.null(norm.factors)) {
                    .self$norm.factors <<- rep(1, length = ncol(count))
                } else {
                    if (length(norm.factors) != ncol(count))
                        stop("\nTCC::ERROR: The length of norm.factors has to be equal to the columns of cuont data.\n")
                    .self$norm.factors <<- norm.factors
                }
                names(norm.factors) <<- norm.factors
            }
            ## Set private argument.
            private$estimated <<- FALSE
            private$simulation <<- FALSE
            private$normalized <<- FALSE
            private$available$norm.method <<- c("tmm", "deseq")
            private$available$test.method <<- data.frame(
                    TwoGroup_NonRep = c(T, T, F, F, F),
                    TwoGroup        = c(T, T, T, T, T),
                    TwoGroup_Paired = c(F, F, F, F, F),
                    MultiGroup      = c(T, T, T, T, T),
                    MultiFactor     = c(T, T, F, T, F),
                    row.names = c("bayseq", "deseq", "ebseq",
                                  "edger", "samseq")
                   )
        }
    )
)

