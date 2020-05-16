module markov_localization


STATIC_ONE_OVER_SQRT_2PI = 1/sqrt(2*pi)
function normpdf(x, mu, std)
    (STATIC_ONE_OVER_SQRT_2PI/std)*exp(-0.5*((x-mu)/std)^2)
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

function initialize_priors(map_size, landmark_positions, position_stdev)
    # set all priors to 0.0
    priors = zeros(map_size)

    # set each landmark positon +/-1 to 1.0/9.0 (9 possible postions)
    norm_term = length(landmark_positions) * (position_stdev * 2 + 1)

    for i=1:length(landmark_positions)
        for j=1:Int(position_stdev)
            priors[(j+landmark_positions[i]+map_size)%map_size] += 1.0/norm_term
            priors[(-j+landmark_positions[i]+map_size)%map_size] += 1.0/norm_term
        end
        priors[landmark_positions[i]] += 1.0/norm_term
    end

  return priors
end



# set standard deviation of control:
control_stdev = 1.0
# set standard deviation of position:
position_stdev = 1.0
# meters vehicle moves per time step
movement_per_timestep = 1.0
# number of x positions on map
map_size = 25

# initialize landmarks
landmark_positions = [6, 11, 21]

# initialize priors
priors = initialize_priors(map_size, landmark_positions, position_stdev)
println(priors)

# step through each pseudo position x (i)
for i = 1:map_size
    pseudo_position = i
    # get the motion model probability for each x position
    motion_prob = motion_model(pseudo_position, movement_per_timestep,
                               priors, map_size, control_stdev);

    print("$pseudo_position $motion_prob\n")
end








end # module
