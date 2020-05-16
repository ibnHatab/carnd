# module markov_localization


STATIC_ONE_OVER_SQRT_2PI = 1/sqrt(2*pi)
function normpdf(x, mu, std)
    (STATIC_ONE_OVER_SQRT_2PI/std)*exp(-0.5*((x-mu)/std)^2)
end

# function to normalize a vector:
function normalize_vector!(a::Array{T,1}) where T
    a  .= a / sum(a)
end


# implement the motion model: calculates prob of being at
# an estimated position at time t
function motion_model(pseudo_position, movement, priors,  map_size, control_stdev)
    # initialize probability
    position_prob = 0.0
    # loop over state space for all possible positions x (convolution):
    for j=1:map_size
        next_pseudo_position = j
        # distance from i to j
        distance_ij = pseudo_position-next_pseudo_position
        # transition probabilities:
        transition_prob = normpdf(distance_ij, movement, control_stdev)
        # estimate probability for the motion model, this is our prior
        position_prob += transition_prob*priors[j]
    end

    return position_prob
end


macro wrap(arr, idx, size) return :($arr[($idx-1+($size))%($size)+1]) end
function initialize_priors(map_size, landmark_positions, position_stdev)
    # set all priors to 0.0
    priors = zeros(map_size)
    # set each landmark positon +/-1 to 1.0/9.0 (9 possible postions)
    norm_term = length(landmark_positions) * (position_stdev * 2 + 1)
    for i=1:length(landmark_positions)
        for j=1:Int(position_stdev)
            priors[(j+landmark_positions[i]-1+(map_size))%(map_size)+1] += 1.0/norm_term
            priors[(-j+landmark_positions[i]-1+(map_size))%(map_size)+1] += 1.0/norm_term
        end
        priors[landmark_positions[i]] += 1.0/norm_term
    end
    return priors
end


function pseudo_range_estimator(landmark_positions, pseudo_position)
    # define pseudo observation vector
    pseudo_ranges = [];
    # loop over number of landmarks and estimate pseudo ranges
    for l=1:length(landmark_positions)
        # estimate pseudo range for each single landmark
        # and the current state position pose_i:
        range_l = landmark_positions[l] - pseudo_position
        # check if distances are positive:
        if (range_l > 0.0)
            push!(pseudo_ranges, range_l)
        end
    end
    # sort pseudo range vector
    sort!(pseudo_ranges)
    return pseudo_ranges;
end


# observation model: calculate likelihood prob term based on landmark proximity
function observation_model(landmark_positions, observations, ranges,
                           distance_max, observation_stdev)
    # initialize observation probability
    distance_prob = 1.0
    reverse!(ranges)
    # run over current observation vector
    for z=1:length(observations)
        # define min distance
        pseudo_range_min = Inf
        # check, if distance vector exists
        if (! isempty(ranges))
            # set min distance
            # remove this entry from ranges-vector
            pseudo_range_min = pop!(ranges)
        else
            # no or negative distances: set min distance to a large number
            pseudo_range_min = Inf
        end
        # estimate the probability for observation model, this is our likelihood
        distance_prob *= normpdf(observations[z], pseudo_range_min, observation_stdev)
    end
    return distance_prob;
end


function main()
    # set standard deviation of control:
    control_stdev = 1.0
    # set standard deviation of position:
    position_stdev = 1.0
    # meters vehicle moves per time step
    movement_per_timestep = 1.0
    # set observation standard deviation
    observation_stdev = 1.0
    # number of x positions on map
    map_size = 25
    # set distance max
    distance_max = map_size
    # initialize landmarks
    landmark_positions = [4, 10, 15, 24]
    # define observations vector, each inner vector represents a set
    #   of observations for a time step
    sensor_obs = [[1,7,12,21], [0,6,11,20], [5,10,19],
                  [4,9,18], [3,8,17], [2,7,16], [1,6,15],
                  [0,5,14], [4,13], [3,12], [2,11], [1,10],
                  [0,9], [8], [7], [6], [5], [4], [3], [2],
                  [1], [0], [], [], []]

    # initialize priors
    priors = initialize_priors(map_size, landmark_positions, position_stdev)
    # initialize posteriors
    posteriors = zeros(map_size)
    # specify time steps
    time_steps = length(sensor_obs)
    # declare observations vector
    observations = Float64[]

    # cycle through time steps
    for t = 1:time_steps
        println("t = $t")
        println("-----Motion----------OBS----------------PRODUCT--")

        if (!isempty(sensor_obs))
            observations = sensor_obs[t]
        else
            observations = Float64(distance_max)
        end

        # step through each pseudo position x(i)
        for i = 1:map_size
            pseudo_position = Float64(i)
            # get the motion model probability for each x position
            motion_prob = motion_model(pseudo_position, movement_per_timestep,
                                       priors, map_size, control_stdev)
            # get pseudo ranges
            ranges = pseudo_range_estimator(landmark_positions, pseudo_position)
            # get observation probability
            observation_prob = observation_model(landmark_positions, observations,
                                                 ranges, distance_max,
                                                 observation_stdev)
            # calculate the ith posterior
            posteriors[i] = motion_prob * observation_prob
            println("$motion_prob \t $observation_prob \t $(motion_prob * observation_prob)")
        end

        println("----------RAW---------------")
        for p = 1:length(posteriors)
            println(posteriors[p])
        end

        # normalize
        normalize_vector!(posteriors)
        println("$(posteriors[t]) \t $(priors[t])")

        println("----------NORMALIZED---------------")
        # update
        priors = posteriors

        for p = 1:length(posteriors)
            println(posteriors[p])
        end
    end
end


main()

#end # module
