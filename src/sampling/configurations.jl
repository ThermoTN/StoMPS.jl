function _hubbard_branch(p,c)
    return c.symmetry == :su2su2 ?
        su2su2_resampling_branch(p.state,su2su2_born_space();
            probability_tolerance=c.probability_tolerance) :
        u1su2_resampling_branch(p.state;probability_tolerance=c.probability_tolerance)
end
_hubbard_rdm(b,i,c) = c.symmetry == :su2su2 ?
    su2su2_resampling_local_rdm(b,i) : u1su2_resampling_local_rdm(b,i)
function _hubbard_project!(b,i,k,c)
    c.symmetry == :su2su2 ?
        su2su2_resampling_project!(b,i,k;probability_tolerance=c.probability_tolerance) :
        u1su2_resampling_project!(b,i,k;probability_tolerance=c.probability_tolerance)
end

"""Evaluate a complete sector configuration under a current normalized MPS."""
function configuration_log_probability(p, configuration,c)
    c.symmetry == :u1su2 && sum(configuration .- 1) != p.N && return -Inf
    b = _hubbard_branch(p,c)
    logp = 0.0
    for (i,k) in enumerate(configuration)
        probability = _hubbard_rdm(b,i,c).sector_probabilities[k]
        probability == 0 && return -Inf
        logp += log(probability)
        _hubbard_project!(b,i,k,c)
    end
    return logp
end

"""Draw an auxiliary configuration by sampling whole local symmetry multiplets."""
function draw_sample_configuration(p,c; rng::AbstractRNG)
    b = _hubbard_branch(p,c)
    X = Int[]
    for i in 1:c.L*c.W
        probabilities = _hubbard_rdm(b,i,c).sector_probabilities
        k = categorical_draw(log.(probabilities);rng)
        push!(X,k); _hubbard_project!(b,i,k,c)
    end
    return X
end
