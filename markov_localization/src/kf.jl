
using LinearAlgebra


function kalman_filter(measurements, x, P)
    # F - state transition matrix
    # H - measurement function
    # R - measurement uncertainty
    for n = 1:length(measurements)
        # measurement update
        Z = [measurements[n]]  # measurement
        y = Z - (H*x)           # motion vector
        S = H*P*H' + R
        K = P*H'*inv(S)
        x = x +(K*y)
        P = (I-(K*H))*P         # uncertanty covariance
        # prediction
        x = (F*x)+u             # estimate
        P = F*P*F'   #
    end
    return x,P
end

############################################
### use the code below to test your filter!
############################################

measurements = [1, 2, 3]

x = [0.0; 0.0] # initial state (location and velocity)
P = [1000.0 0.0;
     0.0    1000.0] # initial uncertainty
u = [0.0; 0.0]  # external motion
F = [1.0 1.0;
     0.0 1.0]   # next state function
H = [1.0 0.0]   # measurement function
R = [1.0]       # measurement uncertainty
I = [1.0 0.0;
     0.0 1.0]   # identity matrix

print(kalman_filter(measurements, x, P))
# output should be:
# x: [[3.9996664447958645], [0.9999998335552873]]
# P: [[2.3318904241194827, 0.9991676099921091], [0.9991676099921067, 0.49950058263974184]]
