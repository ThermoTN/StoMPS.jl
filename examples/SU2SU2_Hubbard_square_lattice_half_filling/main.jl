using StoMPS, MPI

function main(args=ARGS)
    path = isempty(args) ? joinpath(@__DIR__,"config.jl") :
        length(args)==2 && args[1]=="--config" ? abspath(args[2]) :
        error("usage: main.jl [--config PATH]")
    include(path)
    c=Base.invokelatest(config)
    MPI.Init()
    try
        run_application(c;comm=MPI.COMM_WORLD)
    catch exception
        showerror(stderr,exception,catch_backtrace());println(stderr)
        MPI.Comm_size(MPI.COMM_WORLD)>1 && MPI.Abort(MPI.COMM_WORLD,1)
        rethrow()
    finally
        MPI.Finalized() || MPI.Finalize()
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
