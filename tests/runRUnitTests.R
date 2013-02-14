library(RUnit)
library(spatialprobit)

# Define Test Suite for package "tmvtnorm"
testsuite.spatialprobit <- defineTestSuite("spatialprobit", dirs=".", testFileRegexp = "^test.+\\.R")

# Run Test suite
testResult <- runTestSuite(testsuite.spatialprobit)

# print test protocol as text to console or file
printTextProtocol(testResult)
printHTMLProtocol(testResult, fileName = "spatialprobit-unit-test-log.html")