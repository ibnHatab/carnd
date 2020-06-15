module tracking

include("kalman_filter.jl")

@enum SensorType LASER=0 RADAR=1
struct MeasurementPackage
    sensor_type::SensorType
    raw_measurements::Array{Float64, 1}
    timestamp::Int64
end

mutable struct TrackingSt
    # step state
    is_initialized::Bool
    previous_timestamp::Int64
    # acceleration noise components
    noise_ax::Float64
    noise_ay::Float64
    kf::kalman_filter.KalmanFilter
    # ctor
    function TrackingSt()
        # state vector
        x = zeros(4)
        # state covariance matrix
        P = [ 1 0    0    0;
              0 1    0    0;
              0 0 1000    0;
              0 0    0 1000; ]
        # state transistion matrix
        F = [ 1 0 1 0;
              0 1 0 1;
              0 0 1 0;
              0 0 0 1; ]
        # process covariance matrix
        Q = zeros(4,4)
        # measurement matrix
        H = [ 1 0 0 0;
              0 1 0 0; ]
        # measurement covariance matrix
        R = [ 0.0225 0.0;
              0.0    0.0225; ]
        kf = kalman_filter.KalmanFilter(x, P, F, Q, H, R)
        #
        new(false, 0, 5.0, 5.0, kf)
    end
end

function ProcessMeasurement(s::TrackingSt, measurement_pack::MeasurementPackage)
    # Process a single measurement
    if (!s.is_initialized)
        # set the state with the initial location and zero velocity
        s.kf.x = [ measurement_pack.raw_measurements[1]
                   measurement_pack.raw_measurements[2]
                   0
                   0 ]
        s.previous_timestamp = measurement_pack.timestamp
        s.is_initialized = true;
        return
    end

    # compute the time elapsed between the current and previous measurements
    # dt - expressed in seconds
    dt = (measurement_pack.timestamp - s.previous_timestamp) / 1000000.0
    s.previous_timestamp = measurement_pack.timestamp

    dt_2 = dt * dt;
    dt_3 = dt_2 * dt;
    dt_4 = dt_3 * dt;

    # Modify the F matrix so that the time is integrated
    s.kf.F[1, 3] = dt
    s.kf.F[2, 4] = dt

    # set the process covariance matrix Q
    s.kf.Q = [ dt_4/4*s.noise_ax               0  dt_3/2*s.noise_ax               0;
                             0 dt_4/4*s.noise_ay              0   dt_3/2*s.noise_ay;
               dt_3/2*s.noise_ax               0  dt_2*s.noise_ax                 0;
                             0 dt_3/2*s.noise_ay              0     dt_2*s.noise_ay; ]

    #predict

    kalman_filter.Predict(s.kf)

    # measurement update

    kalman_filter.Update(s.kf, measurement_pack.raw_measurements)

    println("x = $(s.kf.x)")
    println("P = $(s.kf.P)")
end

end
