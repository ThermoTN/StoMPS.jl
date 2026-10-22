using FiniteLattices

function polarization_positions(lattice)
    x=[coordinate(lattice,i)[1] for i in 1:size(lattice)]
    return x.-(minimum(x)+maximum(x))/2
end
function polarization_mpo(lattice,positions)
    tree=InteractionTree(size(lattice))
    for i in eachindex(positions)
        addIntr!(tree,U1SU2Fermion.n,i,positions[i];name=:n)
    end
    return AutomataMPO(tree)
end
current_bonds(lattice,positions;t) = [(i,j,positions[j]-positions[i],t) for (i,j) in neighbor(lattice) if i<j]
function current_mpo(lattice,bonds)
    tree=InteractionTree(size(lattice))
    for (i,j,dx,t) in bonds
        addIntr!(tree,U1SU2Fermion.FdagF,(i,j),(true,true),-im*dx*t;
            Z=U1SU2Fermion.Z,name=(:Fdag,:F))
        addIntr!(tree,U1SU2Fermion.FdagF,(j,i),(true,true),im*dx*conj(t);
            Z=U1SU2Fermion.Z,name=(:Fdag,:F))
    end
    return AutomataMPO(tree)
end
function current_observable_tree(lattice,bonds)
    tree=ObservableTree(size(lattice))
    for (i,j,_,_) in bonds
        for (a,b,prefix) in ((i,j,"fwd"),(j,i,"rev"))
            addObs!(tree,U1SU2Fermion.FdagF,(a,b),(true,true);Z=U1SU2Fermion.Z,
                name=(:Fdag,:F),IntrName="$(prefix)_$(i)_$(j)")
        end
    end
    merge!(tree)
    return tree
end
