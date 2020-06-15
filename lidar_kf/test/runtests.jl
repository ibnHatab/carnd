using lidar_kf, Test


function process_measurements(file)
    open(file, "r") do f
        filecontent = read(f,String)
        print(filecontent)
    end
    true
end

@test process_measurements("obj_pose-laser-radar-synthetic-input.txt")
