function make_hamil(vf::Float64,vt::Float64,ωc::Float64,g::Float64,Nx::Int64,Ny::Int64,nmax::Int64)
    kxs = range(-π,π,Nx)
    kys = range(-π,π,Ny)
    ndim = 2 * Nx * Ny * (nmax+1) # 全ヒルベルト空間の次元
    nsect = 2 * Nx * Ny # 光子の個数を確定させたヒルベルト空間の次元
    hamil = spzeros(ComplexF64,ndim,ndim)
    for n in 0:(nmax-1) # n=nmax以外でのハミルトニアンの行列要素を埋める。n=0から始まっていることに注意
        for ix in eachindex(kxs)
            for iy in eachindex(kys)
                kx = kxs[ix]
                ky = kys[iy]
                nk = ix + Nx*(iy-1) # 波数空間のナンバリング,1~Nx*Ny
                hamil[(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2),(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2)] = [ωc*n+vt*kx vf*(kx-1im*ky);
                                                                                                              vf*(kx+1im*ky) ωc*n+vt*kx]
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
            hamil[(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2),(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2)] = [ωc*n+vt*kx vf*(kx-1im*ky);
                                                                                                          vf*(kx+1im*ky) ωc*n+vt*kx]
        end
    end
    return hamil
end

function init_state_ele(Nx::Int64,Ny::Int64)
    kxs = range(-π,π,Nx)
    kys = range(-π,π,Ny)
    nsect = 2 * Nx * Ny
    state_ele = zeros(ComplexF64,nsect) 
    # kx = ky = 0となる波数ラベリングの選定
    ix0 = div(Nx,2) + 1 
    iy0 = div(Ny,2) + 1
    if kxs[ix0] !== 0.0 || kys[iy0] !== 0.0
        error("choose the odd number of Nx and Ny.")
    end
    nk0 = ix0 + Nx*(iy0-1)
    # 初期状態はA sublatticeとB sublatticeに等しく重みづけられた状態になっているとする。
    state_ele[(2*(nk0-1)+1):(2*(nk0-1)+2)] = [1.0, 0.0]
    return state_ele
end

function init_cat(α::ComplexF64,H::SparseMatrixCSC{ComplexF64},Nx::Int64,Ny::Int64,nmax::Int64)
    kxs = range(-π,π,Nx)
    kys = range(-π,π,Ny)
    ndim = 2 * Nx * Ny * (nmax+1)
    nsect = 2 * Nx * Ny
    state = zeros(ComplexF64,ndim)
    state_ele = init_state_ele(Nx,Ny)
    coeff(α,n) = exp(-abs(α)^2/2)*α^n*(1+(-1)^n)/sqrt(2*(1+exp(-2*abs(α)^2))*factorial(n))
    for n in 0:nmax
        state[nsect*n+1:nsect*n+nsect] = coeff(α,n) * state_ele
    end
    return state
end

function init_coherent_p(α::ComplexF64,H::SparseMatrixCSC{ComplexF64},Nx::Int64,Ny::Int64,nmax::Int64)
    kxs = range(-π,π,Nx)
    kys = range(-π,π,Ny)
    ndim = 2 * Nx * Ny * (nmax+1)
    nsect = 2 * Nx * Ny
    state = zeros(ComplexF64,ndim)
    state_ele = init_state_ele(Nx,Ny)
    coeff(α,n) = exp(-abs(α)^2/2)*α^n/sqrt(factorial(n))
    for n in 0:nmax
        state[nsect*n+1:nsect*n+nsect] = coeff(α,n) * state_ele
    end
    return state
end

function init_coherent_m(α::ComplexF64,H::SparseMatrixCSC{ComplexF64},Nx::Int64,Ny::Int64,nmax::Int64)
    kxs = range(-π,π,Nx)
    kys = range(-π,π,Ny)
    ndim = 2 * Nx * Ny * (nmax+1)
    nsect = 2 * Nx * Ny
    state = zeros(ComplexF64,ndim)
    state_ele = init_state_ele(Nx,Ny)
    coeff(α,n) = exp(-abs(α)^2/2)*(-α)^n/sqrt(factorial(n))
    for n in 0:nmax
        state[nsect*n+1:nsect*n+nsect] = coeff(α,n) * state_ele
    end
    return state
end

function make_Jx(vf::Float64,Nx::Int64,Ny::Int64,nmax::Int64)
    kxs = range(-π,π,Nx)
    kys = range(-π,π,Ny)
    ndim = 2 * Nx * Ny * (nmax+1) # 全ヒルベルト空間の次元
    nsect = 2 * Nx * Ny # 光子の個数を確定させたヒルベルト空間の次元
    Jx = spzeros(ComplexF64,ndim,ndim)
    for n in 0:nmax
        for ix in eachindex(kxs)
            for iy in eachindex(kys)
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
    kys = range(-π,π,Ny)
    ndim = 2 * Nx * Ny * (nmax+1) # 全ヒルベルト空間の次元
    nsect = 2 * Nx * Ny # 光子の個数を確定させたヒルベルト空間の次元
    Jy = spzeros(ComplexF64,ndim,ndim)
    for n in 0:nmax
        for ix in eachindex(kxs)
            for iy in eachindex(kys)
                nk = ix + Nx*(iy-1) # 波数空間のナンバリング,1~Nx*Ny
                Jy[(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2),(nsect*n+2*(nk-1)+1):(nsect*n+2*(nk-1)+2)] = [0.0 -1im*vf;
                                                                                                           1im*vf 0.0]
            end
        end
    end
    return Jy
end

function RK4(ψ0::Vector{ComplexF64},H::SparseMatrixCSC{ComplexF64},Δt)
    ψ = deepcopy(ψ0)
    k1 = -1im * Δt * H * ψ
    k2 = -1im * Δt * H * (ψ + k1/2)
    k3 = -1im * Δt * H * (ψ + k2/2)
    k4 = -1im * Δt * H * (ψ + k3)
    ψ += (k1 + 2*k2 + 2*k3 + k4)/6
    normalize!(ψ)
    return ψ
end
