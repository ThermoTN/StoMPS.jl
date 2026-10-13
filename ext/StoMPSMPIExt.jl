module StoMPSMPIExt
using StoMPS, MPI
const SM=StoMPS

struct MPIExecution <: SM.AbstractExecution
    comm::MPI.Comm
end
function SM._execution(comm::MPI.Comm)
    MPI.Initialized() && !MPI.Finalized() || error("initialize MPI before running StoMPS")
    MPIExecution(comm)
end
SM._rank(e::MPIExecution)=MPI.Comm_rank(e.comm)
SM._nprocs(e::MPIExecution)=MPI.Comm_size(e.comm)
SM._gather(x,e::MPIExecution)=MPI.gather(x,e.comm;root=0)
SM._bcast(x,e::MPIExecution)=MPI.bcast(x,e.comm;root=0)
SM._sum(x,e::MPIExecution)=MPI.Allreduce(x,+,e.comm)

# Every rank completes this error collective before the next data collective.
function SM._phase(f,e::MPIExecution)
    result,failure=nothing,""
    try
        result=f()
    catch exception
        failure=sprint(showerror,exception,catch_backtrace())
    end
    failures=MPI.gather(failure,e.comm;root=0)
    message=SM._rank(e)==0 ? join(["rank $(i-1): $s" for (i,s) in enumerate(failures) if !isempty(s)],"\n") : nothing
    message=MPI.bcast(message,e.comm;root=0)
    isempty(message) || error(message)
    return result
end
end
