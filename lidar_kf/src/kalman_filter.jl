module kalman_filter

using LinearAlgebra

mutable struct KalmanFilter
    x::Array{Float64,1}         # state vector
    P::Array{Float64,2}         # state covariance matrix
    F::Array{Float64,2}         # state transistion matrix
    Q::Array{Float64,2}         # process covariance matrix
    H::Array{Float64,2}         # measurement matrix
    R::Array{Float64,2}         # measurement covariance matrix
end

function Predict(kf::KalmanFilter)
    kf.x = kf.F * kf.x
    kf.P = kf.F * kf.P * kf.F' + kf.Q
end

function Update(kf::KalmanFilter, z::AbstractVector)
    z̃ = kf.H * kf.x
    y = z - z̃
    S = kf.H * kf.P * kf.H' + kf.R
    K = kf.P * kf.H' * inv(S)
    # new esimate
    kf.x = kf.x + K*y
    n_ = length(kf.x)
    ID = Matrix{Float64}(I, n_, n_)
    kf.P = (ID - K*kf.H) * kf.P
end

end
