include("workflow_helpers.jl")

@testset "Shared workflow, checkpoint policy and DMRG" begin
    mktempdir() do directory
        c=workflow_config(directory;dynamics_betas=[1/4],dynamics_dt=.05,dynamics_tmax=.1,
            window_cutoff=.1,frequency_list=Float64[])
        result=run_application(c)
        check_workflow(result,c)
        sources=joinpath(directory,"dynamics","sources")
        @test readdir(sources)==["beta_3"]
        dynamics=load_result(joinpath(directory,"results","dynamics_result.jld2"))
        @test dynamics["metadata"]["complete"] && !contains_tensor(dynamics)
        @test all(isfinite,dynamics["observables"]["beta_3"].Cgreater_JP)
        independent=run_dynamics(merge(c,(;dynamics_tmax=.05,window_cutoff=.05,overwrite_output=true)))
        @test independent["observables"]["beta_3"].Cgreater_JP ≈ dynamics["observables"]["beta_3"].Cgreater_JP[1:2]
        cp_path=joinpath(directory,"results","checkpoint.jld2")
        bytes=read(cp_path)
        cp=run_purification(c)
        @test norm(cp.state) ≈ 1
        @test cp.provenance==result["metadata"]["checkpoint_provenance"]
        @test read(cp_path)==bytes
        restarted=run_stochastic(merge(c,(;overwrite_output=true));checkpoint=cp)
        compare_curves(result["observables"]["curve"],restarted["observables"]["curve"])
        ground=run_ground_state(c)
        # Fixed N=2 Hubbard dimer singlet (physical H, no -mu*N).
        @test ground["observables"]["ground_state"].energy_per_site ≈ (c.U-sqrt(c.U^2+16c.t^2))/4 atol=1e-8
        @test load_result(joinpath(directory,"results","ground_state.jld2"))["state"] isa MPS
        @test !contains_tensor(ground)
    end
end

@testset "Full purification reference against Hubbard ED" begin
    mktempdir() do directory
        c=workflow_config(directory;application="SU2SU2_Hubbard_square_lattice_half_filling",
            beta_list=[1/32768,1/8,1/4,1/2],beta0=1/2,environment_disk=false,
            task=:purification_reference,settn_tolerance=1e-12,purification_cutoff=1e-12,
            thermal_seed=20260830,thermal_seed_schedule=:seed_then_seed_plus_beta_index,
            settn_gc_step=false,gc_step=true,gc_sweep=true)
        result=run_purification(c;reference=true,
            before_step=i->Random.seed!(c.thermal_seed+(i==1 ? 0 : i)))
        energies=eigvals(Hermitian(dense_two_site_hubbard(interaction=c.U)))
        rows=result["observables"]["purification"]
        @test result["metadata"]["complete"] && !contains_tensor(result)
        for row in rows
            weights=exp.(-row.beta.*energies)
            @test row.energy_per_site ≈ sum(weights.*energies)/sum(weights)/2 atol=2e-8
            @test row.logZ ≈ log(sum(weights)) atol=2e-8
        end
        cp=load_checkpoint(joinpath(directory,"results","checkpoint.jld2"),c)
        @test cp.logZ0 ≈ last(rows).logZ
        @test last(rows).cumulative_evolution_seconds ≈ sum(r.elapsed_seconds for r in rows)
        @test all(diff([r.wall_elapsed_seconds for r in rows]).>0)
    end
end

@testset "Bulk ratio covariance" begin
    raw=[(;spin_sum=.18,number=.4,double=0.),(;spin_sum=.42,number=.6,double=0.)]
    result=StoMPS.bulk_statistics(raw,zeros(2),ones(Int,2),1,2;probability_tolerance=1e-10)
    @test result.css_tot ≈ 4*.3/(9*.5^2)
    @test result.css_tot_stderr<1e-14
    expanded=StoMPS.bulk_statistics(raw[[1,1,2]],zeros(3),ones(Int,3),1,3;probability_tolerance=1e-10)
    compressed=StoMPS.bulk_statistics(raw,zeros(2),[2,1],1,3;probability_tolerance=1e-10)
    @test collect(values(expanded)) ≈ collect(values(compressed))
end
