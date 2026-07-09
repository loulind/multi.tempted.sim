# want: generate 3 components under ideal scenario, multiply, add noise, and sum
# subject: split into two groups defined by different feature characteristics
# features: 3 modalities, 100 each (can subset later)
# time: S shape, Increasing, Decreasing. Sample equadistant aligned timepoints (10?)
# modality scaling: 
# noise: N(0,1)

library(multi.tempted)

M <- 3
n <- 40
p <- 100
r <- 3
t_min <- 0
t_max <- 1
fn1_s <- function(t){sin(2*pi*t)} # full cycle of sine wave
fn2_incr <- function(t){t-1} # increasing linear function -1 to 1
fn3_decr <- function(t){1-t} # decreasing linear function 1 to -1

