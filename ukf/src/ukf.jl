module ukf

using LinearAlgebra

# input:
# x(k|k) - posterior state
# P(k|k) - posterior covariance
# outpu:
# Xk|k - 2n+1 sigma points
function GenerateSigmaPoints(x::Array{Float64}, P::Array{Float64, 2})
    nₓ = length(x)
    # define spreading parameter
    λ = 3 - nₓ
    # create sigma point matrix
    Xsig = Matrix{Float64}(undef, nₓ, 2*nₓ+1)
    # calculate square root of P
    A = Matrix(cholesky(P).L)
    # set first column of sigma point matrix
    Xsig[:,1] = x
    # set remaining sigma points
    for i = 1:nₓ
        Xsig[:, i+1] = x + sqrt(λ+nₓ) * A[:, i]
        Xsig[:, i+1+nₓ] = x - sqrt(λ+nₓ) * A[:, i]
    end
    Xsig
end

# std_a Process noise standard deviation longitudinal acceleration in m/s^2
# std_yawdd Process noise standard deviation yaw acceleration in rad/s^2
function AugmentedSigmaPoints(x::Array{Float64}, P::Array{Float64, 2};
                              std_a=0.2, std_yawdd=0.2)
    # set state dimension
    nₓ = 5
    # set augmented dimension
    n_aug = 7
    # define spreading parameter
    λ = 3 - n_aug
    # create augmented mean state
    x_aug = vcat(x, [0.0, 0.0])
    # process covariance
    Q = [ std_a*std_a                   0;
                    0 std_yawdd*std_yawdd; ]
    # create augmented covariance matrix
    P_aug = [ P           zeros(5,2);
              zeros(2,5)           Q; ]
    # create sigma point matrix
    Xsig_aug = Matrix(undef, n_aug, 2 * n_aug + 1);
    # create square root matrix
    L = Matrix(cholesky(P_aug).L)
    # create augmented sigma points
    Xsig_aug[:, 1] = x_aug
    for i = 1:n_aug
        Xsig_aug[:, i+1] = x_aug + sqrt(λ+n_aug) * L[:, i]
        Xsig_aug[:, i+1+n_aug] = x_aug - sqrt(λ+n_aug) * L[:, i]
    end
    Xsig_aug
end

function SigmaPointPrediction(Xsig_aug::Array{Float64, 2})
    # set state dimension
    nₓ = 5
    # set augmented dimension
    n_aug = 7
    # create matrix with predicted sigma points as columns
    Xsig_pred = Matrix(undef, nₓ, 2 * n_aug + 1)
    # time diff in sec
    δₜ = 0.1
    # model step function
    function step(x::Array{Float64})
        # predict sigma points
        px, py, v, yaw, yawd, ν_a, ν_yawdd = x
        # predicted state values
        px_p, py_p = 0.0, 0.0
        if (abs(yawd) > 0.001)
            px_p = px + v/yawd * (sin(yaw + yawd*δₜ) - sin(yaw))
            py_p = py + v/yawd * (cos(yaw) - cos(yaw+yawd*δₜ))
        else
            px_p = px + v*δₜ*cos(yaw);
            py_p = py + v*δₜ*sin(yaw);
        end
        v_p = v
        yaw_p = yaw + yawd*δₜ
        yawd_p = yawd
        # add noise
        px_p = px_p + 0.5*ν_a*δₜ*δₜ * cos(yaw)
        py_p = py_p + 0.5*ν_a*δₜ*δₜ * sin(yaw)
        v_p = v_p + ν_a*δₜ
        yaw_p = yaw_p + 0.5*ν_yawdd*δₜ*δₜ
        yawd_p = yawd_p + ν_yawdd*δₜ
        # write predicted sigma point
        [px_p, py_p, v_p, yaw_p, yawd_p]
    end
    # process sigma points
    for i = 1:2*n_aug+1
        x = Xsig_aug[:,i]
        x̃ = step(x)
        Xsig_pred[:, i] = x̃
    end
    Xsig_pred
end

# notmalize radians
function norm_rad(angle)
    while (angle>  pi) angle-=2.0*pi end
    while (angle< -pi) angle+=2.0*pi end
    angle
end

function PredictMeanAndCovariance(Xsig_pred::Array{Float64, 2})
    # set state dimension
    nₓ = 5
    # set augmented dimension
    n_aug = 7
    # define spreading parameter
    λ = 3 - n_aug
    # create vector for weights
    weights = Vector{Float64}(undef, 2*n_aug+1)
    # 2n+1 weights
    weights[1] = λ/(λ+n_aug)
    weights[2:2*n_aug+1] .= 0.5/(n_aug+λ)
    # create vector for predicted state
    x = zeros(nₓ)
    # predicted state mean
    #  iterate over sigma points
    for i = 1:2*n_aug + 1
        x = x + (weights[i] * Xsig_pred[:,i])
    end
    # create covariance matrix for prediction
    P = zeros(nₓ, nₓ)
    for i = 1:2*n_aug + 1
        # state difference
        x_diff = Xsig_pred[:,i] - x;
        x_diff[4] = norm_rad(x_diff[4])
        # angle normalization
        P += weights[i] * (x_diff * x_diff')
    end
    (x, P)
end

function PredictRadarMeasurement(Xsig_pred::Array{Float64, 2})
    # set state dimension
    nₓ = 5
    # set augmented dimension
    n_aug = 7
    # define spreading parameter
    λ = 3 - n_aug
    # set measurement dimension, radar can measure r, phi, and r_dot
    n_z = 3
    # radar measurement noise standard deviation radius in m
    std_radr = 0.3
    # radar measurement noise standard deviation angle in rad
    std_radphi = 0.0175
    # radar measurement noise standard deviation radius change in m/s
    std_radrd = 0.1
    # mean predicted measurement
    z_pred = zeros(n_z)
    # measurement covariance matrix S
    S = zeros(n_z, n_z)
    # create matrix for sigma points in measurement space
    Zsig = Matrix(undef, n_z, 2 * n_aug + 1)
    #
    # set vector for weights
    weights = Vector{Float64}(undef, 2*n_aug+1)
    weights[1] = λ/(λ+n_aug)
    weights[2:2*n_aug+1] .= 0.5/(n_aug+λ)
    # transform sigma points into measurement space
    # 2n+1 simga points
    for i = 1:2 * n_aug + 1
        # extract values for better readability
        p_x = Xsig_pred[1,i]
        p_y = Xsig_pred[2,i]
        v   = Xsig_pred[3,i]
        yaw = Xsig_pred[4,i]
        #
        v1 = cos(yaw)*v
        v2 = sin(yaw)*v
        # measurement model
        Zsig[1,i] = sqrt(p_x*p_x + p_y*p_y)                       # r
        Zsig[2,i] = atan(p_y,p_x)                                 # phi
        Zsig[3,i] = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y)   # r_dot
    end
    # mean predicted measurement
    for i=1:2*n_aug+1
        z_pred = z_pred + weights[i] * Zsig[:,i]
    end
    # innovation covariance matrix S
    # 2n+1 simga points
    for i = 1:2 * n_aug + 1
        # residual
        z_diff = Zsig[:,i] - z_pred
        # angle normalization
        z_diff[2] = norm_rad(z_diff[2])
        #
        S = S + weights[i] * z_diff * z_diff'
    end
    # add measurement noise covariance matrix
    R = [ std_radr*std_radr                      0                   0;
                          0  std_radphi*std_radphi                   0;
                          0                      0 std_radrd*std_radrd;
          ]
    S = S + R
    #
    z_pred, S
end

function UpdateState(Xsig_pred::Array{Float64, 2},
                     x::Array{Float64}, P::Array{Float64, 2},
                     Zsig::Array{Float64, 2}, z_pred::Array{Float64},
                     S::Array{Float64, 2},
                     z::Array{Float64}
                     )
    # set state dimension
    n_x = 5
    # set augmented dimension
    n_aug = 7
    # define spreading parameter
    λ = 3 - n_aug
    # set measurement dimension, radar can measure r, phi, and r_dot
    n_z = 3
    # create vector for weights
    weights = Vector{Float64}(undef, 2*n_aug+1)
    # 2n+1 weights
    weights[1] = λ/(λ+n_aug)
    weights[2:2*n_aug+1] .= 0.5/(n_aug+λ)
    # create matrix for cross correlation Tc
    Tc = zeros(n_x, n_z)
    # calculate cross correlation matrix
    # for 2n+1 simga points
    for i = 1:2 * n_aug + 1
        # residual
        z_diff = Zsig[:,i] - z_pred;
        # angle normalization
        z_diff[2] = norm_rad(z_diff[2])
        # state difference
        x_diff = Xsig_pred[:,i] - x
        # angle normalization
        x_diff[4] = norm_rad(x_diff[4])
        #
        Tc = Tc + weights[i] * x_diff * z_diff'
    end
    # Kalman gain K;
    K = Tc * inv(S)
    # residual
    z_diff = z - z_pred
    # angle normalization
    z_diff[2] = norm_rad(z_diff[2])
    # update state mean and covariance matrix
    x = x + K * z_diff
    P = P - K*S*K'
    #
    x, P
end



end # module
