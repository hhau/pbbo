library(testthat)
library(futile.logger)
library(pbbo)

flog.threshold(DEBUG, name = "pbbo")
flog.threshold(DEBUG)

test_check("pbbo")
