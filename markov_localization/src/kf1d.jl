
# Write a program that will iteratively update and
# predict based on the location measurements
# and inferred motions shown below.

function update(mean1, var1, mean2, var2)
    new_mean = (var2 * mean1 + var1 * mean2) / (var1 + var2)
    new_var = 1.0/(1.0/var1 + 1.0/var2)
    new_mean, new_var
end

function predict(mean1, var1, mean2, var2)
    mean1 + mean2, var1 + var2
end

measurements = [5., 6., 7., 9., 10.]
motion = [1., 1., 2., 1., 1.]
measurement_sig = 4.
motion_sig = 2.
mu = 0.
sig = 10000.

# Insert code here
function kf(mu, sig)
    for (ms, mv) in zip(measurements, motion)
        new_mean, new_var = update(mu, sig, ms, measurement_sig)
        println("update  $new_mean,\t$new_var")
        mu, sig = predict(new_mean, new_var, mv, motion_sig)
        println("predict $mu,\t$sig")
    end
end

kf(mu, sig)
