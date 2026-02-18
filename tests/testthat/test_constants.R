# === Tests for package constants ===================================


# --- TASSEL_JVM constants -----------------------------------------

test_that("TASSEL_JVM is a named list", {
    expect_type(TASSEL_JVM, "list")
    expect_true(length(TASSEL_JVM) > 0)
    expect_true(all(nzchar(names(TASSEL_JVM))))
})

test_that("TASSEL_JVM contains required Java class references", {
    required_keys <- c(
        "ARRAY_LIST",
        "CHROMOSOME",
        "GENOTYPE_TABLE_BUILDER",
        "LOGGING_UTILS",
        "PHENO_BUILDER",
        "R_METHODS",
        "TAXA_LIST_BUILDER"
    )

    for (key in required_keys) {
        expect_true(
            key %in% names(TASSEL_JVM),
            info = paste0("Missing TASSEL_JVM key: ", key)
        )
    }
})

test_that("TASSEL_JVM values are valid Java class paths", {
    for (nm in names(TASSEL_JVM)) {
        val <- TASSEL_JVM[[nm]]
        expect_type(val, "character")
        expect_match(
            val,
            "^[a-zA-Z][a-zA-Z0-9.\\$]*$",
            info = paste0("Invalid Java class path for: ", nm)
        )
    }
})


# --- ANSI constants -----------------------------------------------

test_that("ANSI is a named list with expected keys", {
    expect_type(ANSI, "list")

    expected_keys <- c("BOLD_ON", "BOLD_OFF", "INFO")
    for (key in expected_keys) {
        expect_true(
            key %in% names(ANSI),
            info = paste0("Missing ANSI key: ", key)
        )
    }
})

test_that("ANSI values are character strings", {
    for (nm in names(ANSI)) {
        expect_type(ANSI[[nm]], "character")
        expect_length(ANSI[[nm]], 1)
    }
})

test_that("ANSI escape codes have correct format", {
    expect_match(ANSI$BOLD_ON, "^\033\\[")
    expect_match(ANSI$BOLD_OFF, "^\033\\[")
})

test_that("ANSI INFO is the information symbol", {
    expect_equal(ANSI$INFO, intToUtf8(0x2139))
})


# --- TASSEL_MAVEN constants ---------------------------------------

test_that("TASSEL_MAVEN is a named list", {
    expect_type(TASSEL_MAVEN, "list")
    expect_true(length(TASSEL_MAVEN) > 0)
    expect_true(all(nzchar(names(TASSEL_MAVEN))))
})

test_that("TASSEL_MAVEN contains all required keys", {
    required_keys <- c(
        "BASE_URL",
        "GROUP_PATH",
        "ARTIFACT_ID",
        "VERSION",
        "CLASSIFIER",
        "SHA1_CHECKSUM"
    )

    for (key in required_keys) {
        expect_true(
            key %in% names(TASSEL_MAVEN),
            info = paste0("Missing TASSEL_MAVEN key: ", key)
        )
    }
})

test_that("TASSEL_MAVEN values are non-empty character strings", {
    for (nm in names(TASSEL_MAVEN)) {
        val <- TASSEL_MAVEN[[nm]]
        expect_type(val, "character")
        expect_length(val, 1)
        expect_true(nzchar(val), info = paste0("Empty value for: ", nm))
    }
})

test_that("TASSEL_MAVEN BASE_URL is a valid Maven Central URL", {
    expect_match(TASSEL_MAVEN$BASE_URL, "^https://")
    expect_match(TASSEL_MAVEN$BASE_URL, "maven")
})

test_that("TASSEL_MAVEN VERSION follows semantic versioning", {
    expect_match(TASSEL_MAVEN$VERSION, "^[0-9]+\\.[0-9]+\\.[0-9]+")
})

test_that("TASSEL_MAVEN SHA1_CHECKSUM is a valid 40-character hex string", {
    expect_match(TASSEL_MAVEN$SHA1_CHECKSUM, "^[0-9a-f]{40}$")
})

test_that("TASSEL_MAVEN ARTIFACT_ID is 'tassel'", {
    expect_equal(TASSEL_MAVEN$ARTIFACT_ID, "tassel")
})

test_that("TASSEL_MAVEN CLASSIFIER is 'jar-with-dependencies'", {
    expect_equal(TASSEL_MAVEN$CLASSIFIER, "jar-with-dependencies")
})


