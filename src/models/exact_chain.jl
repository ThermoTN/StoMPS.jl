module OpticalED

using LinearAlgebra

export block_canonical_optical_ed
export finite_time_conductivity


function spin_configurations(length_chain::Int, number_particles::Int)
    0 <= number_particles <= length_chain || throw(ArgumentError(
        "spin-resolved particle number must lie between zero and L",
    ))
    return UInt64[state for state in UInt64(0):((UInt64(1) << length_chain) - 1)
        if count_ones(state) == number_particles]
end

function spin_blocks(length_chain::Int, target_N::Int)
    lower = max(0, target_N - length_chain)
    upper = min(length_chain, target_N)
    blocks = NamedTuple[]
    for number_up in lower:upper
        number_down = target_N - number_up
        number_up <= number_down || continue
        dimension = binomial(length_chain, number_up) *
            binomial(length_chain, number_down)
        multiplicity = number_up == number_down ? 1 : 2
        push!(blocks, (; number_up, number_down, dimension, multiplicity))
    end
    expected = binomial(2 * length_chain, target_N)
    sum(block.multiplicity * block.dimension for block in blocks) == expected ||
        error("spin blocks do not recover the fixed-N Hilbert dimension")
    return blocks
end

function spin_block_operators(length_chain::Int, number_up::Int,
        number_down::Int; hopping::Real=1.0, interaction::Real=8.0)
    up = spin_configurations(length_chain, number_up)
    down = spin_configurations(length_chain, number_down)
    up_index = Dict(state => index for (index, state) in enumerate(up))
    down_index = Dict(state => index for (index, state) in enumerate(down))
    number_up_configurations = length(up)
    dimension = number_up_configurations * length(down)
    hamiltonian = zeros(Float64, dimension, dimension)
    polarization = zeros(Float64, dimension)
    positions = collect(1:length_chain) .- (length_chain + 1) / 2

    for (down_index_value, down_state) in enumerate(down)
        for (up_index_value, up_state) in enumerate(up)
            column = up_index_value +
                (down_index_value - 1) * number_up_configurations
            hamiltonian[column, column] =
                interaction * count_ones(up_state & down_state)
            polarization[column] = sum(positions[site] *
                (((up_state >> (site - 1)) & UInt64(1)) +
                 ((down_state >> (site - 1)) & UInt64(1)))
                for site in 1:length_chain)

            for site in 1:(length_chain - 1)
                bond_mask = (UInt64(1) << (site - 1)) |
                    (UInt64(1) << site)
                left_mask = UInt64(1) << (site - 1)
                right_mask = UInt64(1) << site
                if iszero(up_state & left_mask) != iszero(up_state & right_mask)
                    hopped = up_state ⊻ bond_mask
                    row = up_index[hopped] +
                        (down_index_value - 1) * number_up_configurations
                    hamiltonian[row, column] -= hopping
                end
                if iszero(down_state & left_mask) != iszero(down_state & right_mask)
                    hopped = down_state ⊻ bond_mask
                    row = up_index_value +
                        (down_index[hopped] - 1) * number_up_configurations
                    hamiltonian[row, column] -= hopping
                end
            end
        end
    end
    ishermitian(hamiltonian) || error("spin-block Hamiltonian is not Hermitian")
    return (; hamiltonian, polarization, positions)
end

function parzen_window(time::Real, cutoff::Real)
    cutoff > 0 || throw(ArgumentError("cutoff must be positive"))
    x = abs(Float64(time)) / Float64(cutoff)
    x <= 0.5 && return 1 - 6x^2 + 6x^3
    x <= 1 && return 2(1 - x)^3
    return 0.0
end

function trapezoidal_integral(x, y)
    length(x) == length(y) || throw(DimensionMismatch("x/y length mismatch"))
    length(x) >= 2 || throw(ArgumentError("at least two points are required"))
    return sum((x[index + 1] - x[index]) *
        (y[index] + y[index + 1]) / 2 for index in 1:(length(x) - 1))
end

function finite_time_conductivity(time, correlation, length_chain,
        frequencies; cutoff, eta)
    isapprox(last(time), cutoff; atol=64eps(Float64), rtol=0) ||
        error("ED time grid must end at the transform cutoff")
    kernel = -2 .* imag.(correlation)
    damping = parzen_window.(time, cutoff) .* exp.(-eta .* time)
    damped = kernel .* damping
    conductivity = ComplexF64[
        trapezoidal_integral(time,
            damped .* exp.(1im .* frequency .* time)) / length_chain
        for frequency in frequencies
    ]
    return (; kernel, damping, conductivity)
end

function _add_scaled_block!(log_scale, partition, energy_numerator,
        correlation_numerator, beta_index, block_log_scale,
        block_partition, block_energy, block_correlation)
    if !isfinite(log_scale[beta_index]) || block_log_scale > log_scale[beta_index]
        factor = isfinite(log_scale[beta_index]) ?
            exp(log_scale[beta_index] - block_log_scale) : 0.0
        partition[beta_index] *= factor
        energy_numerator[beta_index] *= factor
        correlation_numerator[beta_index, :] .*= factor
        log_scale[beta_index] = block_log_scale
        block_factor = 1.0
    else
        block_factor = exp(block_log_scale - log_scale[beta_index])
    end
    partition[beta_index] += block_factor * block_partition
    energy_numerator[beta_index] += block_factor * block_energy
    correlation_numerator[beta_index, :] .+=
        block_factor .* block_correlation
    return nothing
end

function block_canonical_optical_ed(length_chain::Int, target_N::Int;
        hopping::Real=1.0, interaction::Real=8.0,
        beta=Float64[1, 2, 4, 8, 16],
        time=collect(0.0:0.2:10.0),
        frequencies=collect(range(0.0, 15.0; length=1201)),
        cutoff::Real=last(time), eta::Real=0.0, on_block=(args...)->nothing)
    length_chain > 0 || throw(ArgumentError("length must be positive"))
    0 <= target_N <= 2length_chain || throw(ArgumentError(
        "target N must lie between zero and 2L",
    ))
    beta = Float64.(beta)
    time = Float64.(time)
    frequencies = Float64.(frequencies)
    all(value -> value > 0 && isfinite(value), beta) ||
        throw(ArgumentError("beta values must be positive and finite"))
    first(time) == 0 && all(isfinite, time) && all(diff(time) .> 0) ||
        throw(ArgumentError("time must start at zero and increase strictly"))
    all(isfinite, frequencies) || throw(ArgumentError("frequencies must be finite"))
    eta >= 0 || throw(ArgumentError("eta must be nonnegative"))

    blocks = spin_blocks(length_chain, target_N)
    number_beta = length(beta)
    number_time = length(time)
    log_scale = fill(-Inf, number_beta)
    partition = zeros(Float64, number_beta)
    energy_numerator = zeros(Float64, number_beta)
    correlation_numerator = zeros(ComplexF64, number_beta, number_time)
    block_records = NamedTuple[]

    for (block_index, block) in enumerate(blocks)
        on_block(:start,(;block_index,block))
        operators = spin_block_operators(length_chain,
            block.number_up, block.number_down; hopping, interaction)
        elapsed_diagonalization = @elapsed decomposition =
            eigen!(Symmetric(operators.hamiltonian))
        eigenvalues = decomposition.values
        eigenvectors = decomposition.vectors

        weighted_vectors = operators.polarization .* eigenvectors
        polarization_eigenbasis = eigenvectors' * weighted_vectors
        transition = abs2.(polarization_eigenbasis)
        @inbounds for column in axes(transition, 2), row in axes(transition, 1)
            transition[row, column] *= eigenvalues[row] - eigenvalues[column]
        end
        phase_minus = cis.(-eigenvalues .* reshape(time, 1, :))
        transition_phase = transition * phase_minus
        phase_plus = conj.(phase_minus)

        for beta_index in eachindex(beta)
            minimum_energy = minimum(eigenvalues)
            shifted_weight = exp.(-beta[beta_index] .*
                (eigenvalues .- minimum_energy))
            block_partition = block.multiplicity * sum(shifted_weight)
            block_energy = block.multiplicity * dot(shifted_weight, eigenvalues)
            weighted_phase = shifted_weight .* phase_plus
            block_correlation = block.multiplicity * 1im .* vec(sum(
                weighted_phase .* transition_phase; dims=1,
            ))
            _add_scaled_block!(log_scale, partition, energy_numerator,
                correlation_numerator, beta_index,
                -beta[beta_index] * minimum_energy,
                block_partition, block_energy, block_correlation)
        end
        push!(block_records, merge(block, (;
            minimum_energy=minimum(eigenvalues),
            maximum_energy=maximum(eigenvalues),
            eigenvalues=copy(eigenvalues),
            elapsed_diagonalization,
        )))
        on_block(:complete,(;block_index,block,elapsed_diagonalization))
        GC.gc()
    end

    log_partition = log_scale .+ log.(partition)
    energy_total = energy_numerator ./ partition
    correlation = correlation_numerator ./ partition
    analyses = NamedTuple[]
    for beta_index in eachindex(beta)
        transformed = finite_time_conductivity(time,
            vec(correlation[beta_index, :]), length_chain, frequencies;
            cutoff, eta)
        push!(analyses, (;
            beta=beta[beta_index],
            temperature=inv(beta[beta_index]),
            log_partition=log_partition[beta_index],
            energy_total=energy_total[beta_index],
            correlation=vec(correlation[beta_index, :]),
            kernel=transformed.kernel,
            conductivity=transformed.conductivity,
        ))
    end
    hilbert_dimension = binomial(2 * length_chain, target_N)
    return (; length=length_chain, target_N, hilbert_dimension,
        hopping=Float64(hopping), interaction=Float64(interaction),
        beta, time, frequencies, cutoff=Float64(cutoff), eta=Float64(eta),
        blocks=block_records, analyses)
end

end # module OpticalED
