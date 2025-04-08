function make_hamil(vf::Float64,ωc::Float64,g::Float64,Nx::Int64,Ny::Int64,nmax::Int64)
    kxs = range(-π,π,Nx)
    kxs = range(-π,π,Ny)
    ndim = 2 * Nx * Ny * (nmax+1) # 全ヒルベルト空間の次元
    nsect = 2 * Nx * Ny # 光子の個数を確定させたヒルベルト空間の次元
    hamil = spzeros(ComplexF64,ndim,ndim)
    for n in 0:(nmax-1) # n=nmax以外でのハミルトニアンの行列要素を埋める。n=0から始まっていることに注意
        for ix in eachindex(kxs)
            for iy in eachindex(kys)
                kx = kxs[ix]
                ky = kys[iy]
                nk = ix + Nx*(iy-1) # 波数空間のナンバリング,1~Nx*Ny
                hamil[(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2),(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2)] = [ωc*n vf*(kx+1im*ky);
                                                                                                              vf*(kx-1im*ky) ωc*n]
                hamil[(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2),(nsect*(n+1)+2*(nk-1)+1):(nsect*(n+1)+2*(nk-1)+2)] = [0.0 -g*sqrt(n+1);
                                                                                                                      -g*sqrt(n+1) 0.0]
                hamil[(nsect*(n+1)+2*(nk-1)+1):(nsect*(n+1)+2*(nk-1)+2),(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2)] = [0.0 -g*sqrt(n+1);
                                                                                                                      -g*sqrt(n+1) 0.0]
            end
        end
    end
    n = nmax # n=nmaxでのハミルトニアンの行列要素を埋める。
    for ix in eachindex(kxs)
        for iy in eachindex(kys)
            kx = kxs[ix]
            ky = kys[iy]
            nk = ix + Nx*(iy-1) # 波数空間のナンバリング,1~Nx*Ny
            hamil[(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2),(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2)] = [ωc*n vf*(kx+1im*ky);
                                                                                                          vf*(kx-1im*ky) ωc*n]
        end
    end
    return hamil
end

function init_cat(α::ComplexF64,Nx::Int64,Ny::Int64,nmax::Int64)
    ndim = 2 * Nx * Ny * (nmax+1)
    state = zeros(ComplexF64,ndim)
    coeff(α,n) = e^(-abs(α)^2/2)*α^n*(1+(-1)^n)/sqrt(2*factorial(n))
    for n in 0:nmax
        for ix in eachindex(kxs)
            for iy in eachindex(kys)
                kx = kxs[ix]
                ky = kys[iy]
                nk = ix + Nx*(iy-1) # 波数空間のナンバリング,1~Nx*Ny
                state[(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2)] = coeff(α,n)*[-1,kx-1im*ky/sqrt(kx^2+ky^2)]
            end
        end
    end
    return state
end

function init_coherent_p(α::ComplexF64,Nx::Int64,Ny::Int64,nmax::Int64)
    ndim = 2 * Nx * Ny * (nmax+1)
    state = zeros(ComplexF64,ndim)
    coeff(α,n) = e^(-abs(α)^2/2)*α^n/sqrt(factorial(n))
    for n in 0:nmax
        for ix in eachindex(kxs)
            for iy in eachindex(kys)
                kx = kxs[ix]
                ky = kys[iy]
                nk = ix + Nx*(iy-1) # 波数空間のナンバリング,1~Nx*Ny
                state[(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2)] = coeff(α,n)*[-1,kx-1im*ky/sqrt(kx^2+ky^2)]
            end
        end
    end
    return state
end

function init_coherent_m(α::ComplexF64,Nx::Int64,Ny::Int64,nmax::Int64)
    ndim = 2 * Nx * Ny * (nmax+1)
    state = zeros(ComplexF64,ndim)
    coeff(α,n) = e^(-abs(α)^2/2)*(-α)^n/sqrt(factorial(n))
    for n in 0:nmax
        for ix in eachindex(kxs)
            for iy in eachindex(kys)
                kx = kxs[ix]
                ky = kys[iy]
                nk = ix + Nx*(iy-1) # 波数空間のナンバリング,1~Nx*Ny
                state[(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2)] = coeff(α,n)*[-1,kx-1im*ky/sqrt(kx^2+ky^2)]
            end
        end
    end
    return state
end

function make_Jx(vf::Float64,Nx::Int64,Ny::Int64,nmax::Int64)
    kxs = range(-π,π,Nx)
    kxs = range(-π,π,Ny)
    ndim = 2 * Nx * Ny * (nmax+1) # 全ヒルベルト空間の次元
    nsect = 2 * Nx * Ny # 光子の個数を確定させたヒルベルト空間の次元
    Jx = spzeros(ComplexF64,ndim,ndim)
    for n in 0:nmax
        for ix in eachindex(kxs)
            for iy in eachindex(kys)
                kx = kxs[ix]
                ky = kys[iy]
                nk = ix + Nx*(iy-1) # 波数空間のナンバリング,1~Nx*Ny
                Jx[(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2),(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2)] = [0.0 -vf;
                                                                                                           -vf 0.0]
            end
        end
    end
    return Jx
end

function make_Jy(vf::Float64,Nx::Int64,Ny::Int64,nmax::Int64)
    kxs = range(-π,π,Nx)
    kxs = range(-π,π,Ny)
    ndim = 2 * Nx * Ny * (nmax+1) # 全ヒルベルト空間の次元
    nsect = 2 * Nx * Ny # 光子の個数を確定させたヒルベルト空間の次元
    Jy = spzeros(ComplexF64,ndim,ndim)
    for n in 0:nmax
        for ix in eachindex(kxs)
            for iy in eachindex(kys)
                kx = kxs[ix]
                ky = kys[iy]
                nk = ix + Nx*(iy-1) # 波数空間のナンバリング,1~Nx*Ny
                Jy[(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2),(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2)] = [0.0 -1im*vf;
                                                                                                           1im*vf 0.0]
            end
        end
    end
    return Jy
end

function RK4(ψ0::Vector{ComplexF64},H::SparseMatrixCSC{ComplexF64,2},Δt)
    ψ = deepcopy(ψ0)
    k1 = -1im * Δt * H * ψ
    k2 = -1im * Δt * H * (ψ + k1/2)
    k3 = -1im * Δt * H * (ψ + k2/2)
    k4 = -1im * Δt * H * (ψ + k3)
    ψ += (k1 + 2*k2 + 2*k3 + k4)/6
    normalize!(ψ)
    return ψ
end
