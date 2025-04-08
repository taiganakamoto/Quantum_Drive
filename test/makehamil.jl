using LinearAlgebra
using SparseArrays
using Random
using HDF5
include("function_hamil.jl")

function main()
    vf = 1.0
    vt = 0.0
    ωc = 1.0
    g = 1e-3
    Nx = 21
    Ny = 21
    nmax = 10
    Δt = 0.01
    T = 50.0
    ts = Vector(0.0:Δt:T)
    Nt = length(ts)
    α = 0.5 + 0.0im
    hamil = make_hamil(vf,vt,ωc,g,Nx,Ny,nmax)
    ψ_cat = init_cat(α,hamil,Nx,Ny,nmax)
    ψ_p = init_coherent_p(α,hamil,Nx,Ny,nmax)
    ψ_m = init_coherent_m(α,hamil,Nx,Ny,nmax)
    Jx = make_Jx(vf,Nx,Ny,nmax)
    Jy = make_Jy(vf,Nx,Ny,nmax)
    Jxs_cat = zeros(Float64,Nt)
    Jys_cat = zeros(Float64,Nt)
    Jxs_p = zeros(Float64,Nt)
    Jys_p = zeros(Float64,Nt)
    Jxs_m = zeros(Float64,Nt)
    Jys_m = zeros(Float64,Nt)
    
    # initial state 
    Jxs_cat[1] = ψ_cat' * Jx * ψ_cat
    Jys_cat[1] = ψ_cat' * Jy * ψ_cat
    Jxs_p[1] = ψ_p' * Jx * ψ_p
    Jys_p[1] = ψ_p' * Jy * ψ_p
    Jxs_m[1] = ψ_m' * Jx * ψ_m
    Jys_m[1] = ψ_m' * Jy * ψ_m

    for nt in 2:Nt
        @show nt
        ψ_cat = RK4(ψ_cat,hamil,Δt)
        Jxs_cat[nt] = ψ_cat' * Jx * ψ_cat
        Jys_cat[nt] = ψ_cat' * Jy * ψ_cat
    end
    for nt in 2:Nt
        @show nt
        ψ_p = RK4(ψ_p,hamil,Δt)
        Jxs_p[nt] = ψ_p' * Jx * ψ_p
        Jys_p[nt] = ψ_p' * Jy * ψ_p
    end
    for nt in 2:Nt
        @show nt
        ψ_m = RK4(ψ_m,hamil,Δt)
        Jxs_m[nt] = ψ_m' * Jx * ψ_m
        Jys_m[nt] = ψ_m' * Jy * ψ_m
    end
    f = h5open("./data_test/test_$(vf)vf$(vt)vt$(ωc)ωc$(g)g$(Nx)Nx$(Ny)Ny$(nmax)nmax$(Δt)Δt$(T)T$(α)α.h5","w")
    write(f,"vf",vf)
    write(f,"vt",vt)
    write(f,"ωc",ωc)
    write(f,"g",g)
    write(f,"Nx",Nx)
    write(f,"Ny",Ny)
    write(f,"nmax",nmax)
    write(f,"Δt",Δt)
    write(f,"ts",ts)
    write(f,"α",α)
    write(f,"Jxs_cat",Jxs_cat)
    write(f,"Jys_cat",Jys_cat)
    write(f,"Jxs_p",Jxs_p)
    write(f,"Jys_p",Jys_p)
    write(f,"Jxs_m",Jxs_m)
    write(f,"Jys_m",Jys_m)
end

main()