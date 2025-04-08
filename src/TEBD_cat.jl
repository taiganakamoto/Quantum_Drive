using ITensors
using ITensorMPS
using HDF5
using Random
using LinearAlgebra
ITensors.Strided.set_num_threads(1)
ITensors.BLAS.set_num_threads(1)
ITensors.enable_threaded_blocksparse(true)

function ITensors.op(::OpName"X", ::SiteType"Qudit", d::Int)
    mat = zeros(d, d)
    for k in 1:(d - 1)
        mat[k + 1, k] = sqrt(k)
        mat[k, k + 1] = sqrt(k)
    end
    return mat
end

function hr(kx::Float64, ky::Float64)
    return kx + 0.0 * ky
end

function hi(kx::Int64, ky::Int64)
    return ky + 0.0 * kx
end

function MPS_ground(kxs::Vector{Float64})
    N = length(kxs)
    sites = siteinds(N)
    ψ = MPS(sites)
    for j in 1:N
        s1 = sites[j-1]
        sp = sites[j]
        H = hr(kxs[j],0.0) * op("N",s1) + hi(0,kxs[j]) * op("Id",s1) + hr(kxs[j],0.0) * op("N",sp) + hi(0,kxs[j]) * op("Id",sp)
        if j == 1
            push!(ψ, H, j)
        else
            push!(ψ, H, j, j-1)
        end
    end
    return ψ
end

function TEBD(ik::Int64,g::Float64,ωc::Float64,sites,tau::Float64)
    N = length(sites)-1
    s1 = sites[j-1]
    sp = sites[j]
    # interaction term
    hj = 
    # photon frequency term
    hj += (Ω/(N-1)) * op("Id",s1) * op("N",sp)
end

function main(Nx::Int64,Δt::Float64,Nt::Int64,md::Int64,cut::Float64)
    kxs = range(-π,π,length=Nx)

    ψ = MPS_ground(kxs)
    sites = siteinds(ψ)

    times = zeros(Float64,Nt+1)
    mds = zeros(Int64,Nt+1)

    # time 0
    step = 0
    time = step*tau
    times[1] = time
    Jxs[1] = #somecode
    Jys[1] = #somecode
    mds[1] = maxlinkdim(ψ)

    for it in 1:Nt
        time = it * Δt
        times[1+it] = time
        #2nd order TEBD
        for ik in 2:length(kxs)
            ψ = swapbondsites(ψ,ik-1)
            sites = siteinds(ψ)
            ψ = apply(TEBD(),ψ;cutoff=cut,maxdim=md)
        end
        for ik in length(kxs):-1:2
            sites = siteinds(ψ)
            ψ = apply(TEBD(),ψ;cutoff=cut,maxdim=md)
            ψ = swapbondsites(ψ,ik-1)
        end
        normalize!(ψ)
        Jxs[1+it] = #somecode
        Jys[1+it] = #somecode
        mds[1+it] = maxlinkdim(ψ)
    end

    # file output
    out = h5open("TEBD_cat.h5","w")
end