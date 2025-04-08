function h(W::Float64,k::Vector{Float64})
    τ1 = [0,1/sqrt(3)]
    τ2 = [-1/2,-1/sqrt(3)/2]
    τ3 = [1/2,-1/sqrt(3)/2]
    return -W * exp(1im*dot(τ1,k)) - W * exp(1im*dot(τ2,k)) - W * exp(1im*dot(τ3,k))
end

function g(W::Float64,A::Float64,k::Vector{Float64})
    τ1 = [0,1/sqrt(3)]
    τ2 = [-1/2,-1/sqrt(3)/2]
    τ3 = [1/2,-1/sqrt(3)/2]
    ex = [1.0,0.0]
    return 1im * W * A * dot(ex,τ1) *exp(1im*dot(τ1,k)) + 1im * W * A * dot(ex,τ2) * exp(1im*dot(τ2,k)) + 1im * W * A * dot(ex,τ3) * exp(1im*dot(τ3,k))
end

function diff_h(W::Float64,k::Vector{Float64})
    τ1 = [0,1/sqrt(3)]
    τ2 = [-1/2,-1/sqrt(3)/2]
    τ3 = [1/2,-1/sqrt(3)/2]
    return 1im * W * τ1 * exp(1im*dot(τ1,k)) + 1im * W * τ2 * exp(1im*dot(τ2,k)) + 1im * W * τ3 * exp(1im*dot(τ3,k))
end

function RK4_pp(h::ComplexF64,g::ComplexF64,α::Float64,ωc::Float64,t::Float64,dt::Float64,σ⁺_pp::ComplexF64,σᶻ_pp::ComplexF64)
    # calculation for the σ⁺_pp and σᶻ_pp
    k⁺1 = -1im * (conj(h) + 2*conj(g)*α*cos(ωc*t)) * σᶻ_pp
    kᶻ1 = -2im * (h + 2*g*α*cos(ωc*t)) * σ⁺_pp + 2im * (conj(h) + 2*conj(g)*α*cos(ωc*t)) * conj(σ⁺_pp)
    k⁺2 = -1im * (conj(h) + 2*conj(g)*α*cos(ωc*(t+dt/2))) * (σᶻ_pp + kᶻ1*dt/2)
    kᶻ2 = -2im * (h + 2*g*α*cos(ωc*(t+dt/2))) * (σ⁺_pp + k⁺1*dt/2) + 2im * (conj(h) + 2*conj(g)*α*cos(ωc*(t+dt/2))) * (conj(σ⁺_pp) + conj(k⁺1)*dt/2)
    k⁺3 = -1im * (conj(h) + 2*conj(g)*α*cos(ωc*(t+dt/2))) * (σᶻ_pp + kᶻ2*dt/2)
    kᶻ3 = -2im * (h + 2*g*α*cos(ωc*(t+dt/2))) * (σ⁺_pp + k⁺2*dt/2) + 2im * (conj(h) + 2*conj(g)*α*cos(ωc*(t+dt/2))) * (conj(σ⁺_pp) + conj(k⁺2)*dt/2)
    k⁺4 = -1im * (conj(h) + 2*conj(g)*α*cos(ωc*(t+dt))) * (σᶻ_pp + kᶻ3*dt)
    kᶻ4 = -2im * (h + 2*g*α*cos(ωc*(t+dt))) * (σ⁺_pp + k⁺3*dt) + 2im * (conj(h) + 2*conj(g)*α*cos(ωc*(t+dt))) * (conj(σ⁺_pp) + conj(k⁺3)*dt)
    σ⁺_pp += + dt/6 * (k⁺1 + 2*k⁺2 + 2*k⁺3 + k⁺4)
    σᶻ_pp += + dt/6 * (kᶻ1 + 2*kᶻ2 + 2*kᶻ3 + kᶻ4)
    return σ⁺_pp, σᶻ_pp
end

function RK4_mm(h::ComplexF64,g::ComplexF64,α::Float64,ωc::Float64,t::Float64,dt::Float64,σ⁺_mm::ComplexF64,σᶻ_mm::ComplexF64)
    # calculation for the σ⁺_mm and σᶻ_mm
    k⁺1 = -1im * (conj(h) - 2*conj(g)*α*cos(ωc*t)) * σᶻ_mm
    kᶻ1 = -2im * (h - 2*g*α*cos(ωc*t)) * σ⁺_mm + 2im * (conj(h) - 2*conj(g)*α*cos(ωc*t)) * conj(σ⁺_mm)
    k⁺2 = -1im * (conj(h) - 2*conj(g)*α*cos(ωc*(t+dt/2))) * (σᶻ_mm + kᶻ1*dt/2)
    kᶻ2 = -2im * (h - 2*g*α*cos(ωc*(t+dt/2))) * (σ⁺_mm + k⁺1*dt/2) + 2im * (conj(h) - 2*conj(g)*α*cos(ωc*(t+dt/2))) * (conj(σ⁺_mm) + conj(k⁺1)*dt/2)
    k⁺3 = -1im * (conj(h) - 2*conj(g)*α*cos(ωc*(t+dt/2))) * (σᶻ_mm + kᶻ2*dt/2)
    kᶻ3 = -2im * (h - 2*g*α*cos(ωc*(t+dt/2))) * (σ⁺_mm + k⁺2*dt/2) + 2im * (conj(h) - 2*conj(g)*α*cos(ωc*(t+dt/2))) * (conj(σ⁺_mm) + conj(k⁺2)*dt/2)
    k⁺4 = -1im * (conj(h) - 2*conj(g)*α*cos(ωc*(t+dt))) * (σᶻ_mm + kᶻ3*dt)
    kᶻ4 = -2im * (h - 2*g*α*cos(ωc*(t+dt))) * (σ⁺_mm + k⁺3*dt) + 2im * (conj(h) - 2*conj(g)*α*cos(ωc*(t+dt))) * (conj(σ⁺_mm) + conj(k⁺3)*dt)
    σ⁺_mm += + dt/6 * (k⁺1 + 2*k⁺2 + 2*k⁺3 + k⁺4)
    σᶻ_mm += + dt/6 * (kᶻ1 + 2*kᶻ2 + 2*kᶻ3 + kᶻ4)
    return σ⁺_mm, σᶻ_mm
end

function RK4_mp(h::ComplexF64,g::ComplexF64,α::Float64,ωc::Float64,t::Float64,dt::Float64,σ⁺_mp::ComplexF64,σ⁻_mp::ComplexF64,σᶻ_mp::ComplexF64)
    # calculation for the σ⁺_mp, σ⁻_mp and σᶻ_mp
    k⁺1 = -1im * (conj(h) - 2im*conj(g)*α*sin(ωc*t)) * σᶻ_mp
    k⁻1 = 1im * (h - 2im*g*α*sin(ωc*t)) * σᶻ_mp
    kᶻ1 = -2im * (h - 2im*g*α*sin(ωc*t)) * σ⁺_mp + 2im * (conj(h) - 2im*conj(g)*α*sin(ωc*t)) * σ⁻_mp
    k⁺2 = -1im * (conj(h) - 2im*conj(g)*α*sin(ωc*(t+dt/2))) * (σᶻ_mp + kᶻ1*dt/2)
    k⁻2 = 1im * (h - 2im*g*α*sin(ωc*(t+dt/2))) * (σᶻ_mp + kᶻ1*dt/2)
    kᶻ2 = -2im * (h - 2im*g*α*sin(ωc*(t+dt/2))) * (σ⁺_mp + k⁺1*dt/2) + 2im * (conj(h) - 2im*conj(g)*α*sin(ωc*(t+dt/2))) * (σ⁻_mp + k⁻1*dt/2)
    k⁺3 = -1im * (conj(h) - 2im*conj(g)*α*sin(ωc*(t+dt/2))) * (σᶻ_mp + kᶻ2*dt/2)
    k⁻3 = 1im * (h - 2im*g*α*sin(ωc*(t+dt/2))) * (σᶻ_mp + kᶻ2*dt/2)
    kᶻ3 = -2im * (h - 2im*g*α*sin(ωc*(t+dt/2))) * (σ⁺_mp + k⁺2*dt/2) + 2im * (conj(h) - 2im*conj(g)*α*sin(ωc*(t+dt/2))) * (σ⁻_mp + k⁻2*dt/2)
    k⁺4 = -1im * (conj(h) - 2im*conj(g)*α*sin(ωc*(t+dt))) * (σᶻ_mp + kᶻ3*dt)
    k⁻4 = 1im * (h - 2im*g*α*sin(ωc*(t+dt))) * (σᶻ_mp + kᶻ3*dt)
    kᶻ4 = -2im * (h - 2im*g*α*sin(ωc*(t+dt))) * (σ⁺_mp + k⁺3*dt) + 2im * (conj(h) - 2im*conj(g)*α*sin(ωc*(t+dt))) * (σ⁻_mp + k⁻3*dt)
    σ⁺_mp += + dt/6 * (k⁺1 + 2*k⁺2 + 2*k⁺3 + k⁺4)
    σ⁻_mp += + dt/6 * (k⁻1 + 2*k⁻2 + 2*k⁻3 + k⁻4)
    σᶻ_mp += + dt/6 * (kᶻ1 + 2*kᶻ2 + 2*kᶻ3 + kᶻ4)
    return σ⁺_mp, σ⁻_mp ,σᶻ_mp
end

function main(W::Float64,A::Float64,α::Float64,ωc::Float64,dt::Float64,Nt::Int64,N::Int64)
    # Set up the parameters
    b1 = [2*π,-2*π/sqrt(3)]
    b2 = [0,4*π/sqrt(3)]
    n1 = Vector(-(N-1)/2:1:(N-1)/2)
    n2 = Vector(-(N-1)/2:1:(N-1)/2)
    n1 /= N
    n2 /= N
    # variables
    σ⁺_pp = zeros(ComplexF64,N,N,Nt)
    σᶻ_pp = zeros(ComplexF64,N,N,Nt)
    σ⁺_mm = zeros(ComplexF64,N,N,Nt)
    σᶻ_mm = zeros(ComplexF64,N,N,Nt)
    σ⁺_mp = zeros(ComplexF64,N,N,Nt)
    σ⁻_mp = zeros(ComplexF64,N,N,Nt)
    σᶻ_mp = zeros(ComplexF64,N,N,Nt)
    ts = zeros(Float64,Nt)
    # Initial conditions
    ts[1] = 0.0
    for i in eachindex(n1)
        for j in eachindex(n2)
            k = b1 * n1[i] + b2 * n2[j]
            σ⁺_pp[i,j,1] = -exp(1im*angle(conj(h(W,k))))/2
            σᶻ_pp[i,j,1] = 0.0
            σ⁺_mm[i,j,1] = -exp(1im*angle(conj(h(W,k))))/2
            σᶻ_mm[i,j,1] = 0.0
            σ⁺_mp[i,j,1] = -exp(1im*angle(conj(h(W,k))))/2 * exp(-2*abs(α)^2)
            σ⁻_mp[i,j,1] = -exp(1im*angle(h(W,k)))/2 * exp(-2*abs(α)^2)
            σᶻ_mp[i,j,1] = 0.0
        end
    end
    # Time evolution
    for step in 2:Nt
        time = dt * (step-1)
        ts[step] = time
        for i in eachindex(n1)
            for j in eachindex(n2)
                k = b1 * n1[i] + b2 * n2[j]
                h_val = h(W,k)
                g_val = g(W,A,k)
                σ⁺_pp[i,j,step], σᶻ_pp[i,j,step] = RK4_pp(h_val,g_val,α,ωc,time,dt,σ⁺_pp[i,j,step-1],σᶻ_pp[i,j,step-1])
                σ⁺_mm[i,j,step], σᶻ_mm[i,j,step] = RK4_mm(h_val,g_val,α,ωc,time,dt,σ⁺_mm[i,j,step-1],σᶻ_mm[i,j,step-1])
                σ⁺_mp[i,j,step], σ⁻_mp[i,j,step], σᶻ_mp[i,j,step] = RK4_mp(h_val,g_val,α,ωc,time,dt,σ⁺_mp[i,j,step-1],σ⁻_mp[i,j,step-1],σᶻ_mp[i,j,step-1])
            end
        end
    end
    # Save the results
    f = h5open("../data/Heisenberg_cat.h5", "w")
    write(f, "σ⁺_pp", σ⁺_pp)
    write(f, "σᶻ_pp", σᶻ_pp)
    write(f, "σ⁺_mm", σ⁺_mm)
    write(f, "σᶻ_mm", σᶻ_mm)
    write(f, "σ⁺_mp", σ⁺_mp)
    write(f, "σ⁻_mp", σ⁻_mp)
    write(f, "σᶻ_mp", σᶻ_mp)
    write(f, "ts", ts)
    write(f, "n1", n1)
    write(f, "n2", n2)
    write(f, "b1", b1)
    write(f, "b2", b2)
end