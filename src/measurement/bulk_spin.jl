# Physical bulk selection, scheduled spin measurements and ensemble ratios.
const CORRELATION_DEFINITION = (
    name=:manuscript_bulk_internal_css_tot,
    spin_operator="physical S_i dot S_j; SU2-invariant thermal state",
    pair_rule="unique undirected nearest-neighbor bonds with both endpoints in bulk",
    formula="4*mean(spin_sum)/(9*Nbulk*mean(single_probability)^2)",
    normalization="square of bulk site-averaged thermal single occupancy",
    errorbar="first-order ratio propagation with covariance",
)

function correlation_schedule(c)
    betas = sample_betas(c)
    schedule = NamedTuple[]
    for T in c.correlation_temperatures
        isfinite(T) && T > 0 || error("correlation temperatures must be finite and positive")
        matches = findall(b -> isapprox(b, inv(T); rtol=64eps(Float64), atol=0), betas)
        length(matches)==1 || error("correlation temperature $T requires exactly one saved sample beta=$(inv(T)); edit beta_list explicitly")
        i = only(matches)
        push!(schedule, (;index=i, beta=betas[i], requested_temperature=Float64(T),
            temperature=inv(betas[i]), temperature_over_t=inv(betas[i])/c.t))
    end
    allunique(r.index for r in schedule) || error("duplicate correlation temperatures")
    return sort!(schedule; by=r->r.index)
end

function bulk_geometry(c, model)
    columns = c.bulk_columns
    columns isa AbstractUnitRange{<:Integer} && !isempty(columns) &&
        1 <= first(columns) <= last(columns) <= c.L ||
        error("bulk_columns must be an explicit nonempty integer range within 1:L")
    sites = findall(i -> model.lattice[i][1] in columns, 1:c.L*c.W)
    selected = Set(sites)
    bonds = filter(ij -> ij[1] in selected && ij[2] in selected, model.bonds)
    isempty(bonds) && error("bulk contains no internal nearest-neighbor bonds")
    return (;columns=collect(columns), sites, bonds, Nbulk=length(sites),
        number_bonds=length(bonds), spatial_normalizer=3length(sites))
end

function bulk_operators(c, model, bulk)
    nsites = c.L*c.W
    number_tree, double_tree = InteractionTree(nsites), InteractionTree(nsites)
    for i in bulk.sites
        addIntr!(number_tree, U1SU2Fermion.n, i, 1.0; name=:n)
        addIntr!(double_tree, U1SU2Fermion.nd, i, 1.0; name=:nd)
    end
    return (;spin=hubbard_spin_correlation_mpo(nsites,bulk.bonds),
        number=AutomataMPO(number_tree), double=AutomataMPO(double_tree))
end

function bulk_sample_observables(state, operators; disk=true)
    return (;spin_sum=hubbard_expectation(state,operators.spin;disk),
        number=hubbard_expectation(state,operators.number;disk),
        double=hubbard_expectation(state,operators.double;disk))
end

function bulk_statistics(raw, log_weight, multiplicity, Nbulk, Nv; probability_tolerance)
    spin = [r.spin_sum for r in raw]
    number = [r.number for r in raw]
    double = [r.double for r in raw]
    single = (number .- 2 .* double) ./ Nbulk
    all(x -> -probability_tolerance <= x <= 1+probability_tolerance, single) ||
        error("invalid bulk single-occupancy probability")
    stats(x) = multiplicity_weighted_statistics(x,log_weight,multiplicity,Nv)
    q, n, d, p = stats(spin), stats(number), stats(double), stats(single)
    p.mean > probability_tolerance || error("bulk single occupancy too small for css_tot normalization")
    factor = 4/(9Nbulk)
    css_tot = factor*q.mean/p.mean^2
    # The same sample supplies numerator and denominator; their covariance is
    # retained by propagating the centered joint influence.
    influence = factor .* ((spin .- q.mean)./p.mean^2 .-
        2q.mean .* (single .- p.mean)./p.mean^3)
    css_tot_stderr = stats(influence).errorbar
    isfinite(css_tot) && isfinite(css_tot_stderr) || error("invalid normalized spin correlation")
    return (;css_tot,css_tot_stderr,spin_sum=q.mean,spin_sum_stderr=q.errorbar,
        bulk_number=n.mean,bulk_number_stderr=n.errorbar,
        bulk_double=d.mean,bulk_double_stderr=d.errorbar,
        bulk_single_mean=p.mean,bulk_single_stderr=p.errorbar,
        bulk_density=n.mean/Nbulk,bulk_doping=n.mean/Nbulk-1,
        bulk_doping_stderr=n.errorbar/Nbulk,normalization=p.mean^2,
        Nbulk,spatial_normalizer=3Nbulk)
end

function correlation_measurements(particles,s,c,e,stats,operators,bulk,index,out)
    rows = _phase(e) do
        values=zeros(Float64,length(s),3)
        for p in particles
            v=bulk_sample_observables(p.state,operators;disk=c.environment_disk)
            values[findfirst(q->q.id==p.id,s),:]=[v.spin_sum,v.number,v.double]
        end
        values
    end
    rows = _sum(rows,e)
    raw = [(;id=s[i].id,spin_sum=rows[i,1],number=rows[i,2],double=rows[i,3])
        for i in eachindex(s)]
    # Save raw moments before evaluating the normalized correlation.
    result = _broadcast(e) do
        out.data["diagnostics"]["pending_correlations"]=(;index,raw,scalars=s,bulk)
        result=bulk_statistics(raw,stats.w,stats.m,bulk.Nbulk,stats.Nv;
            probability_tolerance=c.probability_tolerance)
        delete!(out.data["diagnostics"],"pending_correlations")
        result
    end
    return (;result...,raw)
end
