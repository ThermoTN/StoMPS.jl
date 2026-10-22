"""A normalized Born sample with its configuration, observables and sampling weights."""
Base.@kwdef mutable struct HubbardSample{S<:DenseMPS}
    id::Int
    sample_index::Int = id
    lineage::Int
    generation::Int = 0
    state::S
    configuration::Vector{Int}
    N::Int
    energy::Float64
    a::Float64 = 0.0
    h::Float64 = 0.0
    m::Int = 1
    log_p0::Float64
    log_Q::Float64
    compression_norm::Float64
    diagnostics = nothing
end

sample_betas(c) = c.beta_list[c.beta_list .>= c.beta0]

"""Draw/project the explicitly selected Hubbard Born law, then compress."""
function hubbard_born(parent, c; rng::AbstractRNG, configuration=nothing)
    _validate_hubbard_ensemble(c)
    kw = (; trunc=tensor_truncation(c.sample_bond_dimension,c.born_cutoff),
        probability_tolerance=c.probability_tolerance)
    if c.symmetry == :su2su2
        return configuration === nothing ?
            su2su2_purification_born_sample(parent,su2su2_born_space();rng,kw...) :
            su2su2_purification_born_project(parent,su2su2_born_space(),configuration;kw...)
    elseif c.ensemble == :canonical
        return configuration === nothing ?
            u1su2_fixed_n_purification_born_sample(parent,Int(c.Ne);rng,kw...) :
            u1su2_fixed_n_purification_born_project(parent,Int(c.Ne),configuration;kw...)
    elseif c.ensemble == :grand_canonical
        return configuration === nothing ?
            u1su2_purification_born_sample(parent;rng,kw...) :
            u1su2_purification_born_project(parent,configuration;kw...)
    end
end

function hubbard_sample(born, model, c, id; generation=0, lineage=id, m=1,
        log_Q=born.log_probability)
    N = c.symmetry == :su2su2 ? c.L*c.W : sum(born.configuration .- 1)
    c.ensemble == :canonical && N != c.Ne && error("Born particle number mismatch")
    return HubbardSample(;id,lineage,generation,state=born.state,
        configuration=copy(born.configuration),N,
        energy=hubbard_expectation(born.state,model.H;disk=c.environment_disk),
        a=born.log_probability-log_Q,m,log_p0=born.log_probability,log_Q,
        compression_norm=born.compression_norm)
end

"""Evolve one Born sample using an explicitly seeded evolution step."""
function advance_hubbard_sample!(p, model, c, index; step!::Function)
    beta = sample_betas(c)
    diagnostics = step!(p.state,model.H,-(beta[index]-beta[index-1])/2,
        c.sample_bond_dimension,index,p.generation,p.lineage)
    # Use the step's squared-norm log increment, or normalize a raw step here.
    increment=hasproperty(diagnostics,:log_norm_increment) ? diagnostics.log_norm_increment :
        normalize_logweight!(p.state)
    p.a += increment; p.h += increment
    p.energy = hubbard_expectation(p.state,model.H;disk=c.environment_disk)
    p.diagnostics=hasproperty(diagnostics,:maximum_truncation) ? diagnostics : nothing
    return diagnostics
end

sample_scalars(p::HubbardSample) = (;
    p.id,p.sample_index,p.lineage,p.generation,p.configuration,p.N,p.energy,p.a,p.h,p.m,
    p.log_p0,p.log_Q,p.compression_norm)
