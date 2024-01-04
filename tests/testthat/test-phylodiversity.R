

sample <- rsample(1, phylocom$phylo)
test_that('mpd',
{
  expect_equal(mpd(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted = TRUE), expected = c(4.250000, 4.944444, 5.833333, 6.944444, 7.750000, 7.109375)
    , tolerance = 1e-5)
  expect_vector(mpd(sample, cophenetic(phylocom$phylo), abundance.weighted = TRUE))
  expect_true(all(mpd(sample, cophenetic(phylocom$phylo), abundance.weighted = TRUE) > 0))
  expect_no_error(mpd(sample, cophenetic(phylocom$phylo), abundance.weighted = TRUE), message = NULL, class = NULL)
  expect_no_warning(mpd(sample, cophenetic(phylocom$phylo), abundance.weighted = TRUE), message = NULL, class = NULL)
  expect_no_message(mpd(sample, cophenetic(phylocom$phylo), abundance.weighted = TRUE), message = NULL, class = NULL)
})


test_that(
  'mntd', {
  expect_vector(mntd(sample, cophenetic(phylocom$phylo), abundance.weighted = TRUE))
  expect_true(all(mntd(sample, cophenetic(phylocom$phylo), abundance.weighted = TRUE) > 0))
  expect_no_error(mntd(sample, cophenetic(phylocom$phylo), abundance.weighted = TRUE), message = NULL, class = NULL)
  expect_no_warning(mntd(sample, cophenetic(phylocom$phylo), abundance.weighted = TRUE), message = NULL, class = NULL)
  expect_no_message(mntd(sample, cophenetic(phylocom$phylo), abundance.weighted = TRUE), message = NULL, class = NULL)
}
)

test_that(
  'ses.mpd', {
  expect_vector(ses.mpd(phylocom$sample, cophenetic(phylocom$phylo), null.model = "taxa.labels"))
  expect_no_error(ses.mpd(phylocom$sample, cophenetic(phylocom$phylo), null.model = "taxa.labels"), message = NULL, class = NULL)
  expect_no_warning(ses.mpd(phylocom$sample, cophenetic(phylocom$phylo), null.model = "taxa.labels"), message = NULL, class = NULL)
  expect_no_message(ses.mpd(phylocom$sample, cophenetic(phylocom$phylo), null.model = "taxa.labels"), message = NULL, class = NULL)
}
)


test_that(
  'ses.mntd', {
  expect_vector(ses.mntd(phylocom$sample, cophenetic(phylocom$phylo), null.model = "taxa.labels"))
  expect_no_error(ses.mntd(phylocom$sample, cophenetic(phylocom$phylo), null.model = "taxa.labels"), message = NULL, class = NULL)
  expect_no_warning(ses.mntd(phylocom$sample, cophenetic(phylocom$phylo), null.model = "taxa.labels"), message = NULL, class = NULL)
  expect_no_message(ses.mntd(phylocom$sample, cophenetic(phylocom$phylo), null.model = "taxa.labels"), message = NULL, class = NULL)
}
)

test_that(
  'pd', {
  expect_vector(pd(phylocom$sample, phylocom$phylo))
  expect_true(all(pd(phylocom$sample, phylocom$phylo) > 0))
  expect_no_error(pd(phylocom$sample, phylocom$phylo), message = NULL, class = NULL)
  expect_no_warning(pd(phylocom$sample, phylocom$phylo), message = NULL, class = NULL)
  expect_no_message(pd(phylocom$sample, phylocom$phylo), message = NULL, class = NULL)
}
)

test_that(
  'ses.pd', {
  expect_vector(ses.pd(phylocom$sample, phylocom$phylo, null.model = "taxa.labels", runs = 99))
  expect_no_error(ses.pd(phylocom$sample, phylocom$phylo, null.model = "taxa.labels", runs = 99), message = NULL, class = NULL)
  expect_no_warning(ses.pd(phylocom$sample, phylocom$phylo, null.model = "taxa.labels", runs = 99), message = NULL, class = NULL)
  expect_no_message(ses.pd(phylocom$sample, phylocom$phylo, null.model = "taxa.labels", runs = 99), message = NULL, class = NULL)
}
)


