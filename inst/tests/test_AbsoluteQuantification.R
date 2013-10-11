library(testthat)
data(UPS2MS)

# AbsoluteQuantification.default
test_that("AbsoluteQuantification.default", {
	expect_that(AbsoluteQuantification.default(UPS2_SRM)[c("is_calibrated","mfe","r.squared","calibration_covar")],equals(list("is_calibrated"=TRUE,"mfe"=1.670342,"r.squared"=0.9512209,"calibration_covar"=16.9878), tolerance = .001))
	expect_that(AbsoluteQuantification.default(UPS2_LFQ)[c("is_calibrated","mfe","r.squared","calibration_covar")],equals(list("is_calibrated"=TRUE,"mfe"=2.120428,"r.squared"=0.8075852,"calibration_covar"=11.07064), tolerance = .001))
	expect_that(AbsoluteQuantification.default(UPS2_SC)[c("is_calibrated","mfe","r.squared","calibration_covar")],equals(list("is_calibrated"=TRUE,"mfe"=1.974945,"r.squared"=0.8753872,"calibration_covar"=78.7692), tolerance = .001))
})

# predict.AbsoluteQuantification
test_that("predict.AbsoluteQuantification", {
	expect_that(predict.AbsoluteQuantification(AbsoluteQuantification.default(UPS2_SRM))$prediction$concentration[1:5],equals(c(0.29041894,4.77261553,1.57177383,5.84919886,5.46618250), tolerance = .001))
	expect_that(predict.AbsoluteQuantification(AbsoluteQuantification.default(UPS2_LFQ))$prediction$concentration[1:5],equals(c(4.327546,5.516119,5.921030,5.955162,3.742530), tolerance = .001))
	expect_that(predict.AbsoluteQuantification(AbsoluteQuantification.default(UPS2_SC))$prediction$concentration[1:5],equals(c(4.189661,7.047456,5.431220,5.999418,3.261362), tolerance = .001))
})

# folderror.AbsoluteQuantification
test_that("folderror.AbsoluteQuantification", {
	expect_that(folderror.AbsoluteQuantification(c(10,9,8),c(5,4,2)),equals(c(2.00,2.25,4.00)))
})

# cval.AbsoluteQuantification
test_that("cval.AbsoluteQuantification: Monte Carlo", {
	set.seed(131)
	expect_that(cval.AbsoluteQuantification(AbsoluteQuantification.default(UPS2_SRM),cval_method="mc",mcx=2)$cval[c("r.squared","mfe")],equals(list("r.squared"=0.9292566,"mfe"=1.761073), tolerance = .001))

	set.seed(131)
	expect_that(cval.AbsoluteQuantification(AbsoluteQuantification.default(UPS2_LFQ),cval_method="mc",mcx=2)$cval[c("r.squared","mfe")],equals(list("r.squared"=0.5395614,"mfe"=3.251972), tolerance = .001))

	set.seed(131)
	expect_that(cval.AbsoluteQuantification(AbsoluteQuantification.default(UPS2_SC),cval_method="mc",mcx=2)$cval[c("r.squared","mfe")],equals(list("r.squared"=0.6703485,"mfe"=2.294331), tolerance = .001))
})
