using JLD
using LinearAlgebra
using LorentzDrudeMetals
using ProgressMeter
using Random
BLAS.set_num_threads(3)
include("./fdfd/FiniteDifferenceHomogenisation.jl")

sf = 100 # Factor for lamdas (100, 11points), (10, 101points), (1, 1001points)
v_diameter = [10, 15, 20, 25, 30]
gap = 2
v_res = [1.00, 0.75]

for diameter in v_diameter
    for res in v_res
        s_res = 1000*res
        begin
            lam0s = union(0.5:0.001*sf:1.50) # microns
            metal = LorentzDrudeMetals.Au
            dname = joinpath("./data/Au$(diameter)_$(s_res)/")
            @showprogress for lam0 in shuffle(lam0s)
                s_lam0 = 1000*lam0
                k0 = 2pi/lam0/1000 # [1/nm]
                fname = joinpath(dname, "lam0=$(s_lam0).jld")
                if isfile(fname)
                    continue
                end
                ep, mu, dx, dy, dz = makesphereshex([diameter/2,gap/2], [metal[lam0],1], res)
                neffs = getneff(ep, mu, dx, dy, dz, k0)
                save(fname, "neffs", neffs, "lam0", lam0, "dx", dx, "dy", dy, "dz", dz)
            end
        end
    end
end
