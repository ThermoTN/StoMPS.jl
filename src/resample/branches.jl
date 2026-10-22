"""Process-local whole-sector measurement branch."""
mutable struct SectorBranch{S,P,Q}
    state::S
    physical_space::P
    sector_list::Q
    next_site::Int
    active::Bool
    last_local
end

function _sector_branch(snapshot,physical_space,sector_list;probability_tolerance)
    probability_tolerance>0 || throw(ArgumentError("probability tolerance must be positive"))
    isapprox(norm(snapshot),1;atol=probability_tolerance,rtol=0) ||
        error("resampling requires a normalized snapshot")
    all(t->t isa MPSTensor{4} && codomain(t)[2]==physical_space,snapshot) ||
        error("resampling requires rank-4 tensors with the selected physical space")
    state=deepcopy(snapshot)
    Center(state)==[1,1] || canonicalize!(state,1;trunc=notrunc())
    SectorBranch(state,physical_space,sector_list,1,true,nothing)
end

function _sector_local_rdm(branch::SectorBranch,site::Integer)
    branch.active && site==branch.next_site || error("inactive or out-of-order measurement branch")
    center=branch.state[site].A
    @tensor rdm[physical;physical_prime] :=
        center[left physical;marker right]*conj(center[left physical_prime;marker right])
    coefficients=Float64[real(only(block(rdm,s))) for s in branch.sector_list]
    dimensions=Int[dim(_born_sector_space(branch.physical_space,s)) for s in branch.sector_list]
    probabilities=coefficients.*dimensions
    all(isfinite,probabilities) && all(>=(0),probabilities) || error("invalid sector probabilities")
    result=(;rdm,rdm_coefficients=coefficients,sector_probabilities=probabilities,
        trace_residual=abs(real(tr(rdm))-1),hermiticity_residual=norm(rdm-rdm'),
        probability_residual=abs(sum(probabilities)-1),sector_list=branch.sector_list,
        marker_dimensions=dimensions)
    branch.last_local=result
    return result
end

function _sector_project!(branch::SectorBranch,site::Integer,sector_index::Integer;
        probability_tolerance=1e-10)
    branch.active || return (likelihood=0.0,active=false,probability_residual=0.0)
    site==branch.next_site && branch.last_local!==nothing || error("measure this site before projection")
    likelihood=branch.last_local.sector_probabilities[sector_index]
    branch.next_site=site+1
    branch.last_local=nothing
    if iszero(likelihood)
        branch.active=false
        return (;likelihood,active=false,probability_residual=0.0)
    end
    marker=_born_sector_space(branch.physical_space,branch.sector_list[sector_index])
    embedding=isometry(branch.physical_space,marker)
    projector=embedding*embedding'
    center=branch.state[site].A
    @tensor projected[left physical;marker right] :=
        projector[physical;input]*center[left input;marker right]
    branch.state[site]=MPSTensor(projected/sqrt(likelihood))
    site<length(branch.state) && canonicalize!(branch.state,site+1;trunc=notrunc())
    residual=abs(norm(branch.state)-1)
    residual<=probability_tolerance || error("sector projection lost normalization")
    return (;likelihood,active=true,probability_residual=residual)
end

function _validate_su2su2_resampling_space(space::GradedSpace)
    space==su2su2_born_space() || throw(ArgumentError("require canonical su2su2_born_space()"))
end

function su2su2_resampling_branch(snapshot,physical_space::GradedSpace;probability_tolerance=1e-10)
    _validate_su2su2_resampling_space(physical_space)
    _sector_branch(snapshot,physical_space,collect(sectors(physical_space));probability_tolerance)
end
function u1su2_resampling_branch(snapshot;probability_tolerance=1e-10)
    physical_space=codomain(snapshot[1])[2]
    _validate_u1su2_space(physical_space)
    _sector_branch(snapshot,physical_space,_u1su2_sector_list();probability_tolerance)
end
su2su2_resampling_local_rdm(branch,site)=_sector_local_rdm(branch,site)
u1su2_resampling_local_rdm(branch,site)=_sector_local_rdm(branch,site)
su2su2_resampling_project!(branch,site,index;kwargs...)=_sector_project!(branch,site,index;kwargs...)
u1su2_resampling_project!(branch,site,index;kwargs...)=_sector_project!(branch,site,index;kwargs...)
