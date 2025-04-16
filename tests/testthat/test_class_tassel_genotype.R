test_that("TasselGenotype class construction and methods work", {
    # Mock a Java reference object
    mockJObj <- structure(list(), class = "jobjRef")
    
    # Create a TasselGenotype object
    gt <- new("TasselGenotype",
              dispData = list(name = "test_genotype"),
              jRefObj = mockJObj,
              jMemAddress = "0x123",
              jClass = "net.maizegenetics.dna.snp.GenotypeTable")
    
    # Test class structure
    expect_s4_class(gt, "TasselGenotype")
    expect_type(gt@dispData, "list")
    expect_s3_class(gt@jRefObj, "jobjRef")
    expect_type(gt@jMemAddress, "character")
    expect_type(gt@jClass, "character")
    
    # Test javaRefObj method
    expect_identical(javaRefObj(gt), mockJObj)
    
    # Test show method (capturing output)
    expect_output(show(gt), "TasselGenotype")  # Basic check that output contains class name
})

test_that("readGenotype handles invalid inputs correctly", {
    # Test invalid file path
    expect_error(readGenotype("nonexistent/path.txt"))
    
    # Test invalid input type
    expect_error(readGenotype(list()), "Unsupported data type")
})