

test_that('mpd',
{
  data(phylocom)
  data("dummy_normal")
  data("dummy_all_NA.rda")
  data("dummy_one_na.rda")
  
  
  expect_equal(mpd(dummy_normal, cophenetic(phylocom$phylo), abundance.weighted=TRUE),  expected = c(7.954789, 7.895544, 7.978408, 7.946582, 7.917633, 7.987820)
    ,tolerance = 1e-5 )
  ##  expect_identical(mpd(dummy_all_NA, cophenetic(phylocom$phylo), abundance.weighted=TRUE),  expected = c(NA, NA, NA, NA, NA, NA))
  expect_vector(mpd(dummy_normal, cophenetic(phylocom$phylo), abundance.weighted=TRUE))
  expect_equal(mpd(dummy_one_na, cophenetic(phylocom$phylo), abundance.weighted=TRUE),  expected = c(7.864198, 7.666667, 7.695502, 7.975207, 7.697778, 7.857143)
    ,tolerance = 1e-5 )
  expect_true(all(mpd(dummy_normal, cophenetic(phylocom$phylo), abundance.weighted=TRUE)>0))
  expect_no_error(mpd(dummy_normal, cophenetic(phylocom$phylo), abundance.weighted=TRUE), message = NULL, class = NULL)
  expect_no_warning(mpd(dummy_normal, cophenetic(phylocom$phylo), abundance.weighted=TRUE), message = NULL, class = NULL)
  expect_no_message(mpd(dummy_normal, cophenetic(phylocom$phylo), abundance.weighted=TRUE), message = NULL, class = NULL)
})



test_that(
  'mntd', {
  expect_vector(mntd(dummy_normal, cophenetic(phylocom$phylo), abundance.weighted=TRUE))
  expect_true(all(mntd(dummy_normal, cophenetic(phylocom$phylo), abundance.weighted=TRUE)>0))
  expect_no_error(mntd(dummy_normal, cophenetic(phylocom$phylo), abundance.weighted=TRUE), message = NULL, class = NULL)
  expect_no_warning(mntd(dummy_normal, cophenetic(phylocom$phylo), abundance.weighted=TRUE), message = NULL, class = NULL)
  expect_no_message(mntd(dummy_normal, cophenetic(phylocom$phylo), abundance.weighted=TRUE), message = NULL, class = NULL)
}
)

test_that(
  'ses.mpd', {
  expect_vector(ses.mpd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels") )
 ## expect_true(all(ses.mpd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels")>0))
  expect_no_error(ses.mpd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels"), message = NULL, class = NULL)
  expect_no_warning(ses.mpd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels"), message = NULL, class = NULL)
  expect_no_message(ses.mpd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels"), message = NULL, class = NULL)
}
)

test_that(
'ses.mntd', {
expect_vector(ses.mntd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels") )
## expect_true(all( ses.mntd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels")>0))
expect_no_error( ses.mntd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels"), message = NULL, class = NULL)
expect_no_warning(ses.mntd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels"), message = NULL, class = NULL)
expect_no_message(ses.mntd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels"), message = NULL, class = NULL)
}
)

test_that(
'pd', {
expect_vector(pd(phylocom$sample, phylocom$phylo))
expect_true(all(pd(phylocom$sample, phylocom$phylo) >0))
expect_no_error(pd(phylocom$sample, phylocom$phylo) , message = NULL, class = NULL)
expect_no_warning(pd(phylocom$sample, phylocom$phylo), message = NULL, class = NULL)
expect_no_message(pd(phylocom$sample, phylocom$phylo), message = NULL, class = NULL)
}
)

test_that(
'ses.pd', {
expect_vector(ses.pd(phylocom$sample, phylocom$phylo, null.model="taxa.labels", runs=99))
##expect_true(all(ses.pd(phylocom$sample, phylocom$phylo, null.model="taxa.labels", runs=99)>0))
expect_no_error(ses.pd(phylocom$sample, phylocom$phylo, null.model="taxa.labels", runs=99) , message = NULL, class = NULL)
expect_no_warning(ses.pd(phylocom$sample, phylocom$phylo, null.model="taxa.labels", runs=99), message = NULL, class = NULL)
expect_no_message(ses.pd(phylocom$sample, phylocom$phylo, null.model="taxa.labels", runs=99), message = NULL, class = NULL)
}
)


