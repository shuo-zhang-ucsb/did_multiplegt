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

did_multiplegt(df, Y, G, T, D, controls, placebo = placebo, dynamic = dynamic,
               brep = 2, cluster = "nr", covariance = TRUE, average_effect = "prop_number_switchers") #add placebo and dynamics

cluster = NULL
recat_treatment = NULL
trends_nonparam = NULL
trends_lin = "nr"
did_multiplegt(df, Y, G, T, D, controls, placebo = placebo, dynamic = dynamic,
               trends_lin = trends_lin, trends_nonparam = trends_nonparam) #add placebo and dynamics


cluster = NULL
recat_treatment = NULL
trends_nonparam = "nr"
trends_lin = NULL
did_multiplegt(df, Y, G, T, D, controls, placebo = placebo, dynamic = dynamic,
               trends_lin = trends_lin, trends_nonparam = trends_nonparam) #add placebo and dynamics
