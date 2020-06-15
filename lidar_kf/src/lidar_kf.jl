module lidar_kf

using DelimitedFiles

include("tracking.jl")

function process_measurements()
    # Set Measurements
    in_file_name = "../test/obj_pose-laser-radar-synthetic-input.txt";
    measurement_pack_list = []
    # scan file
    mix_array = readdlm(in_file_name, '\t')
    for i = 1:3
        p = mix_array[i,:]
        type, x, y, timestamp = p[1:4]
        if type == "L"
            meas_package = tracking.MeasurementPackage(tracking.LASER, [x,y], timestamp)
            push!(measurement_pack_list, meas_package)
        end
    end
    # Create a Tracking instance
    st = tracking.TrackingSt()
    # call the ProcessingMeasurement() function for each measurement
    # start filtering from the second frame
    # (the speed is unknown in the first frame)
    println(measurement_pack_list)
    for k = 1:length(measurement_pack_list)
        tracking.ProcessMeasurement(st, measurement_pack_list[k])
    end
end

process_measurements()

end # module
