test_that("enrichSTRING_API retries on transient HTTP errors and eventually fails", {
    # Track how many times POST is called
    call_count <- 0L

    # Mock httr::POST to always return a 503 response
    mock_response <- structure(list(), class = "response")

    local_mocked_bindings(
        POST = function(...) {
            call_count <<- call_count + 1L
            mock_response
        },
        status_code = function(resp) 503L,
        content = function(resp, ...) "Service Unavailable",
        .package = "httr"
    )

    expect_error(
        enrichSTRING_API(
            genes       = c("TP53", "BRCA1"),
            max_retries = 2L,
            retry_delay = 0
        ),
        "STRING API request failed"
    )

    # Should have tried initial attempt + 2 retries = 3 total calls
    expect_equal(call_count, 3L)
})

test_that("enrichSTRING_API does not retry on non-retryable errors", {
    call_count <- 0L
    mock_response <- structure(list(), class = "response")

    local_mocked_bindings(
        POST = function(...) {
            call_count <<- call_count + 1L
            mock_response
        },
        status_code = function(resp) 400L,
        content = function(resp, ...) "Bad Request",
        .package = "httr"
    )

    expect_error(
        enrichSTRING_API(
            genes       = c("TP53", "BRCA1"),
            max_retries = 3L,
            retry_delay = 0
        ),
        "STRING API request failed"
    )

    # Should only call POST once for a non-retryable error
    expect_equal(call_count, 1L)
})

test_that("enrichSTRING_API succeeds after initial transient failures", {
    call_count <- 0L
    mock_response <- structure(list(), class = "response")

    tsv_text <- paste(
        "term\tdescription\tp_value\tfdr\tcategory\tnumber_of_genes\tnumber_of_genes_in_background\tncbiTaxonId\tinputGenes",
        "hsa04110\tCell cycle\t0.001\t0.01\tKEGG\t5\t100\t9606\tTP53",
        sep = "\n"
    )

    local_mocked_bindings(
        POST = function(...) {
            call_count <<- call_count + 1L
            mock_response
        },
        status_code = function(resp) if (call_count < 3L) 503L else 200L,
        content = function(resp, as = "text", ...) {
            if (call_count < 3L) "Service Unavailable" else tsv_text
        },
        .package = "httr"
    )

    result <- enrichSTRING_API(
        genes           = c("TP53", "BRCA1"),
        category        = "KEGG",
        adjpvalueCutoff = 1.0,
        max_retries     = 3L,
        retry_delay     = 0
    )

    expect_equal(call_count, 3L)
    expect_s4_class(result, "enrichResult")
})
