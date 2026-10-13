@testset "JP dynamics against canonical ED" begin
    c=workflow_config(nothing;dynamics_dt=.05,dynamics_tmax=.1,window_cutoff=.1,
        frequency_list=[0.,1.],broadening=0.,environment_disk=false)
    model=hubbard_model(c)
    cp=run_purification(c)
    trace=zeros(ComplexF64,3)
    for (id,X) in enumerate(([1,3],[2,2],[3,1]))
        born=hubbard_born(cp.state,c;rng=MersenneTwister(1),configuration=X)
        p=hubbard_sample(born,model,c,id)
        response=StoMPS.conductivity_particle(p,c,model,1)
        trace .+= exp(born.log_probability).*response.Cgreater_JP
        @test response.t0_residual_abs<1e-8
    end
    exact=StoMPS.OpticalED.block_canonical_optical_ed(2,2;hopping=c.t,interaction=c.U,
        beta=[c.beta0],time=[0.,.05,.1],frequencies=c.frequency_list,cutoff=.1,eta=0.)
    @test trace ≈ exact.analyses[1].correlation atol=3e-7
end
