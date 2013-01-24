#
# RUnit test cases for function LeSagePaceExperiment()
#
#
test.LeSagePaceExperiment <- function() {
  # 1. Check for valid inputs
  checkException(LeSagePaceExperiment(n=100, beta=1))
  checkException(LeSagePaceExperiment(n=100, beta="A"))
}