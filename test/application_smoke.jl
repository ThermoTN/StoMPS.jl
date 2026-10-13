# Example environment: OUTPUT [APPLICATION_NAME] [BASELINE_OUTPUT].
using Test, MPI, StoMPS
include("workflow_helpers.jl")
application=length(ARGS)>=2 ? ARGS[2] : "U1SU2_Hubbard_square_lattice_canonical_ensemble"
c=workflow_config(abspath(ARGS[1]);application)
MPI.Init()
try
    result=run_application(c;comm=MPI.COMM_WORLD)
    # A failure on one rank must reach every rank before the next collective.
    failed=false
    try
        StoMPS._phase(StoMPS._execution(MPI.COMM_WORLD)) do
            MPI.Comm_rank(MPI.COMM_WORLD)==MPI.Comm_size(MPI.COMM_WORLD)-1 && error("injected phase failure")
        end
    catch
        failed=true
    end
    all_failed=MPI.Allreduce(failed,&,MPI.COMM_WORLD)
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @testset "MPI workflow and numerical output" begin
            @test all_failed
            check_workflow(result,c)
            if length(ARGS)>=3
                baseline=StoMPS.JLD2.load(joinpath(ARGS[3],"history.jld2"))["history"]
                compare_curves(baseline,result["observables"]["curve"])
            end
        end
    end
catch exception
    showerror(stderr,exception,catch_backtrace())
    MPI.Comm_size(MPI.COMM_WORLD)>1 && MPI.Abort(MPI.COMM_WORLD,1)
    rethrow()
finally
    MPI.Finalized() || MPI.Finalize()
end
