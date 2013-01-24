library(RUnit)
library(spatialprobit)

# command line arguments
Args <- commandArgs()

# working directory
wd=Args[5]
if (is.null(wd) | is.na(wd)) wd="."
setwd(wd)

setwd("F:/R/spatialprobit")
library(spatialprobit)

# Define Test Suite for package "tmvtnorm"
testsuite.spatialprobit <- defineTestSuite("spatialprobit", dirs="tests", testFileRegexp = "^test.+\\.R")

# Run Test suite
testResult <- runTestSuite(testsuite.spatialprobit)

# print test protocol as text to console or file
printTextProtocol(testResult)
printHTMLProtocol(testResult, fileName = "spatialprobit-unit-test-log.html")