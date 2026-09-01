  function moving_average(vector::Vector{T}, window::Int) where T<:Number
    n = length(vector)
    out = Vector{Float64}(undef, n - window + 1)
    
    # Calculate the first window sum
    current_sum = sum(@view vector[1:window])
    out[1] = current_sum / window
    
    # Roll through the rest of the array by adding the next element and subtracting the oldest
    for i in 2:(n - window + 1)
        current_sum += vector[i + window - 1] - vector[i - 1]
        out[i] = current_sum / window
    end
    return out
end

  function NoUnloadBuckling_And_Uniform_TimeInc(npoints,raw_λ,raw_P,raw_times)
    @assert length(raw_times) == length(raw_λ) == length(raw_P)
    P_max, max_pos = findmax(raw_P)
    # Remove the near zero/negative stresses after unloading
    last_measure = findlast(σ -> σ > 0.1*P_max, raw_P)
    last_time = raw_times[last_measure]
    # Generate a new time series and interpolate data
    times = range(raw_times[1], last_time, length=npoints)
    interp_λ = LinearInterpolation(raw_λ, raw_times)
    interp_P = LinearInterpolation(raw_P, raw_times)
    λ = interp_λ(times)
    P = interp_P(times)

    # Generate the new struct
    λ_max = round(maximum(λ); digits=1)
    Δt = times[11] - times[10]
    λ_max, Δt, λ, P, times
  end

  function OnlyLoad_OnlyIncreasingDispl_MovingAvg_UniformTimeInc(npoints,window,raw_λ,raw_P,raw_times)
    @assert length(raw_times) == length(raw_λ) == length(raw_P)
    P_max, max_pos = findmax(raw_P)
    raw_λ = raw_λ[1:max_pos]
    raw_P = raw_P[1:max_pos]
    raw_times = raw_times[1:max_pos]
    λ_max = raw_λ[end]
    # remove noisy non-increasing stretch
    increasing_indices = [1; findall(diff(raw_λ) .> 1e-6) .+ 1]
    raw_λ = raw_λ[increasing_indices]
    raw_P = raw_P[increasing_indices]
    raw_times = raw_times[increasing_indices]

    #moving averages
    #raw_λ = moving_average(raw_λ,window)
    #raw_P = moving_average(raw_P,window)
    #raw_times = moving_average(raw_times,window)

    # Generate a new time series and interpolate data
    times = range(raw_times[1], raw_times[end], length=npoints)
    interp_λ = LinearInterpolation(raw_λ, raw_times)
    interp_P = LinearInterpolation(raw_P, raw_times)
    λ = interp_λ(times)
    P = interp_P(times)

    Δt = times[2] - times[1]
    λ_max, Δt, λ, P, times
  end