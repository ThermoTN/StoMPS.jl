abstract type AbstractExecution end
struct SerialExecution <: AbstractExecution end
_execution(::Nothing)=SerialExecution()
_execution(e::AbstractExecution)=e
_execution(comm)=error("load MPI and initialize it before passing a communicator")
_rank(::SerialExecution)=0
_nprocs(::SerialExecution)=1
_phase(f,::SerialExecution)=f()
_gather(x,::SerialExecution)=[x]
_bcast(x,::SerialExecution)=x
_sum(x,::SerialExecution)=x
_root(f,e::AbstractExecution)=_phase(()->_rank(e)==0 ? f() : nothing,e)
_broadcast(f,e::AbstractExecution)=_bcast(_root(f,e),e)
owns(id,e::AbstractExecution)=mod(id-1,_nprocs(e))==_rank(e)
rng(c,domain,event,item)=MersenneTwister(deterministic_seed(c.master_seed,domain,event,item))

function global_scalars(particles,e)
    parts=_gather(sample_scalars.(particles),e)
    _broadcast(()->sort!(reduce(vcat,parts);by=p->p.id),e)
end
