using LinearAlgebra
using SparseArrays
using Random
using HDF5
include("function_hamil.jl")

function main()
    vf = 1.0
    ωc = 1.0
    g = 1e-3
    Nx = 10
    Ny = 10
    nmax = 5
    Δt = 0.01
    Nt = 100
    α = 2.0
    hamil = make_hamil(vf,ωc,g,Nx,Ny,nmax)
    ψ_cat = init_cat(α,Nx,Ny,nmax)
    ψ_p = init_coherent_p(α,Nx,Ny,nmax)
    ψ_m = init_coherent_m(α,Nx,Ny,nmax)
    Jx = make_Jx(vf,Nx,Ny,nmax)
    Jy = make_Jy(vf,Nx,Ny,nmax)
    Jxs_cat = zeros(Float64,Nt)
    Jys_cat = zeros(Float64,Nt)
    Jxs_p = zeros(Float64,Nt)
    Jys_p = zeros(Float64,Nt)
    Jxs_m = zeros(Float64,Nt)
    Jys_m = zeros(Float64,Nt)
    for nt in 1:Nt
        ψ_cat = RK4(ψ_cat,hamil,Δt)
        Jxs_cat[nt] = ψ' * Jx * ψ
        Jys_cat[nt] = ψ' * Jy * ψ
    end
    for nt in 1:Nt
        ψ_p = RK4(ψ_p,hamil,Δt)
        Jxs_p[nt] = ψ' * Jx * ψ
        Jys_p[nt] = ψ' * Jy * ψ
    end
    for nt in 1:Nt
        ψ_m = RK4(ψ_m,hamil,Δt)
        Jxs_m[nt] = ψ' * Jx * ψ
        Jys_m[nt] = ψ' * Jy * ψ
    end
    f = h5open("./test.h5","w")
    write(f,"vf",vf)
    write(f,"ωc",ωc)
    write(f,"g",g)
    write(f,"Nx",Nx)
    write(f,"Ny",Ny)
    write(f,"nmax",nmax)
    write(f,"Δt",Δt)
    write(f,"Nt",Nt)
    write(f,"α",α)
    write(f,"Jxs_cat",Jxs_cat)
    write(f,"Jys_cat",Jys_cat)
    write(f,"Jxs_p",Jxs_p)
    write(f,"Jys_p",Jys_p)
    write(f,"Jxs_m",Jxs_m)
    write(f,"Jys_m",Jys_m)
end

main()