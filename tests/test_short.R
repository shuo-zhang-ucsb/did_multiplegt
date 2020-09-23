library(DIDmultiplegt)
library(wooldridge)

df <- wagepan
Y = "lwage"
G = "nr"
T = "year"
D = "union"
controls = c("hours")
placebo = 2
dynamic = 2

result = did_multiplegt(df, Y, G, T, D, controls, placebo = placebo, dynamic = dynamic,
                        brep = 2, cluster = "nr", covariance = TRUE, average_effect = "prop_number_switchers")

assertthat::are_equal(result$placebo_2, -0.9930194)
assertthat::are_equal(result$N_placebo_2, 2158)
assertthat::are_equal(result$placebo_1, 0.08446127)
assertthat::are_equal(result$N_placebo_1,2842)
assertthat::are_equal(result$effect, 0.02147226)
assertthat::are_equal(result$N_effect, 3815)
assertthat::are_equal(result$N_switchers_effect, 508)


cluster = NULL
recat_treatment = NULL
trends_nonparam = NULL
trends_lin = "nr"
result = did_multiplegt(df, Y, G, T, D, controls, placebo = placebo, dynamic = dynamic,
                        trends_lin = trends_lin, trends_nonparam = trends_nonparam)


cluster = NULL
recat_treatment = NULL
trends_nonparam = "black"
trends_lin = NULL
result = did_multiplegt(df, Y, G, T, D, controls, placebo = placebo, dynamic = dynamic,
                        trends_lin = trends_lin, trends_nonparam = trends_nonparam)
