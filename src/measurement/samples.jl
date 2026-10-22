function sample_logweights(s, c, lambda, lambda0)
    return [p.a + (c.ensemble == :grand_canonical ? (lambda-lambda0)*p.N : 0)
        for p in s]
end

function sample_statistics(s,c,lambda,lambda0)
    w = sample_logweights(s,c,lambda,lambda0)
    m = [p.m for p in s]
    Nv = sum(m)
    metrics = multiplicity_metrics(w,m,Nv)
    eta_p = clamp(metrics.eta_p*(Nv/c.Ns),1/c.Ns,length(s)/c.Ns)
    eta = c.eta_weight == :estimator ? metrics.eta :
        sampling_efficiency([p.h for p in s])
    energy = c.ensemble == :grand_canonical ?
        fitted_statistics([p.energy for p in s],[p.N for p in s],w,m) :
        multiplicity_weighted_statistics([p.energy for p in s],w,m,Nv)
    errors=conditional_weight_errors(w,m)
    return (;w,m,Nv,eta,eta_p,energy=energy.mean,
        lnZ_error=errors.lnZ_error,eta_error=c.eta_weight==:estimator ? errors.eta_error : NaN,
        stderr=energy.errorbar,Np=length(Set(Tuple(p.configuration) for p in s)),
        K=length(s),mean_N=logweighted_mean([p.N for p in s],w .+ log.(m)))
end

"""Estimate log-partition and sampling-efficiency errors from paired moments."""
function conditional_weight_errors(log_weight,multiplicity)
    Nv=sum(multiplicity)
    Nv<=1 && return (;lnZ_error=NaN,eta_error=NaN)
    w=exp.(log_weight.-maximum(log_weight))
    mu1=sum(multiplicity.*w)/Nv
    mu2=sum(multiplicity.*w.^2)/Nv
    sem(f)=sqrt(max(0.0,sum(multiplicity.*abs2.(f))/(Nv*(Nv-1))))
    (;lnZ_error=sem((w.-mu1)./mu1),
        eta_error=sem((2mu1/mu2).*(w.-mu1).-(mu1^2/mu2^2).*(w.^2 .-mu2)))
end
