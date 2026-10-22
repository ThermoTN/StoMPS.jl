"""Expectation with an application-selected disk environment."""
function hubbard_expectation(state, operator; disk=true)
    with_environment(state,operator;disk) do env
        value = scalar!(env; normalize=true)
        isfinite(value) && abs(imag(value)) <= 1e-8 || error("invalid expectation")
        return Float64(real(value))
    end
end

"""Site-resolved number and double occupancy for a U1SU2 state."""
function hubbard_local_observables(state; disk=true)
    number, double = Float64[], Float64[]
    for i in 1:length(state)
        for (op, values) in ((U1SU2Fermion.n,number),(U1SU2Fermion.nd,double))
            tree = InteractionTree(length(state))
            addIntr!(tree, op, i, 1.0)
            push!(values,hubbard_expectation(state,AutomataMPO(tree); disk))
        end
    end
    return (; number, double)
end
