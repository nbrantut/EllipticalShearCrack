using QuadGK

"""
    ellipticalshearcrack(t, c, θ, ϕ, Σ, N, cdcs; kwargs...)

Compute a given component of either normalised acceleration or stress rate around a dynamically expanding, self-similar, elliptical shear crack, using the formulae given in Richards, IJSS 1973.

Normalisations follow those of Richards 1973. Numerical integration is done with the QuadGK package (note: this method seems to be immune to nearby singularities in the integration path, so that the explicit extraction of those singularities as described by Richards in the appendix seem a bit redundant. Here I implemented the full method for safety but tests show it might be overkill. Richards (pers. comm.) used Romberg integration so that singularities were more of a problem.)

# Arguments:
 - `t`: an array of normalised time
 - `c`: a String indicating which component to compute. Values can be ''ux'', ''uy'', ''uz'' (normalised accelerations) or ''zx'', ''zy'', ''zz'' (normalised stress rates).
 - `θ`: angle of point of interest from normal to shear crack.
 - `ϕ`: azimuth of point of interest from x direction.
 - `Σ`: rupture speed in x direction.
 - `N`: rupture speed in y direction.
 - `cdcs`: ratio of P to S wave speed.
 - additional keyword arguments can be passed to the `quadgk` numerical integration algorithm.

# Acknowledgments

I am deeply grateful to Paul Richards for his help and for sharing his original Fortran code.

"""
function ellipticalshearcrack(t,c,θ,ϕ,Σ,N,cdcs;kwargs...)

    p = param(θ,ϕ,Σ,N,cdcs)

    # pseudoarrivals
    G = p[:D]*sin(θ)^2*(p[:F] + p[:D]*cos(θ)^2) - Σ^2*N^2
    if G>0
        w0d = √(p[:D]*cos(θ)^2*(1 - p[:D]*τα('d',p)^2*sin(θ)^2)/G)
        τ0d = √(w0d^2*Σ^2*N^2 + p[:D])/(sin(θ)* p[:D])
        w0s = √(p[:D]*cos(θ)^2*(1 - p[:D]*τα('s',p)^2*sin(θ)^2)/G)
        τ0s = √(w0s^2*Σ^2*N^2 + p[:D])/(sin(θ)* p[:D])
        #println("Pseudoarrivals at τ0d=",round(τ0d,digits=2))
        #println("                  τ0s=",round(τ0s,digits=2))
    else
        τ0d = -1.
        τ0s = -1.
    end
    

    # arrivals
    τd = τα('d',p)
    τs = τα('s',p)

    sol = Float64[]
    #solR = Float64[]
    #solI = Float64[]

    wd = 1+1.0im
    ws = 1+1.0im

    for τ in t
        I = 0.
        R = 0.

        if τ>τd
            Td = √(τ^2 - τd^2)
            if (G>0) && (τ>τ0d)
                wd = wα(τ,'d',wd,p)
                I += quadgk(x -> integrand_smooth_nonsingular(Td, x, wd, τ, 'd', c ,p), 0, π/2;kwargs...)[1] + integralofSα(Td,wd,τ,'d',c,p)
                #I += quadgk(x -> integrand_smooth(Td, x, τ, 'd', c ,p), 0, π/2;kwargs...)[1]
                R += residue('d',c,wd,p)
            else
                I += quadgk(x -> integrand_smooth(Td, x, τ, 'd', c ,p), 0, π/2;kwargs...)[1]
            end
        end
        if τ>τs
            Ts = √(τ^2 - τs^2)
            if (G>0) && (τ>τ0s)
                ws = wα(τ,'s',ws,p)
                I += quadgk(x -> integrand_smooth_nonsingular(Ts, x, ws, τ, 's', c ,p), 0, π/2;kwargs...)[1] + integralofSα(Ts,ws,τ,'s',c,p)
                #I += quadgk(x -> integrand_smooth(Ts, x, τ, 's', c ,p), 0, π/2;kwargs...)[1]
                R += residue('s',c,ws,p)
            else
                I += quadgk(x -> integrand_smooth(Ts, x, τ, 's', c ,p), 0, π/2;kwargs...)[1]
            end
        end
        push!(sol, I+R)
        #push!(solR, R)
        #push!(solI, I)

    end
    return sol
end

function param(θ,ϕ,Σ,N,cdcs)
    p = Dict()
    p[:θ] = θ
    p[:ϕ] = ϕ
    p[:Σ] = Σ
    p[:N] = N
    p[:cdcs] = cdcs
    p[:D] = Σ^2*cos(ϕ)^2 + N^2*sin(ϕ)^2
    p[:F] = Σ^2*sin(ϕ)^2 + N^2*cos(ϕ)^2
    p[:Δ] = abs((Σ^2 - N^2)*cos(ϕ)*sin(ϕ))

    return p
end

function τα(α,p)
    return (α=='d' ? 1.0 : p[:cdcs])
end

function qα(w,t,α,p)
    return √(t^2 - w^2 - τα(α,p)^2)*cos(p[:θ]) + 1.0im*t*sin(p[:θ])
end

function dqαdt(w,t,α,p)
    return t*cos(p[:θ])/√(t^2 + w^2 - τα(α,p)^2) + 1.0im*sin(p[:θ])
end

function dqαdw(w,t,α,p)
    return -w*cos(p[:θ])/√(t^2 + w^2 - τα(α,p)^2)
end

function d2qαdw2(w,t,α,p)
    sq  = √(t^2 + w^2 - τα(α,p)^2)
    return cos(p[:θ])/sq*(-1 -w^2/sq^2)
end

function qσν(w,p)
    return (w*p[:Δ] + 1.0im*√(w^2*p[:Σ]^2*p[:N]^2 + p[:D]))/p[:D]
end

function dqσνdw(w,p)
    return (p[:Δ] + 1.0im*(w*p[:Σ]^2*p[:N]^2)/√(w^2*p[:Σ]^2*p[:N]^2 + p[:D]))/p[:D]
end

function wα(τ,α,w0,p)
    τa = τα(α,p)
    
    Z(w) = -1.0im*qσν(w,p)*sin(p[:θ])+ √(qσν(w,p)^2 + w^2 + τa^2)*cos(p[:θ])
    dZdw(w) = -1.0im*sin(p[:θ])*dqσνdw(w,p) + cos(p[:θ])/√(qσν(w,p)^2 + w^2 + τa^2)*(qσν(w,p)*dqσνdw(w,p) + w)

    # find wa(tau) using newton-raphson
    tol = 1e-8
    imax = 500
    counter = 0
    wa = w0
    residual = τ - Z(wa)
    while abs(residual)>tol
        wa = wa + residual/dZdw(wa)
        residual = τ - Z(wa)
        counter += 1
        if counter > imax
            error("Max. number of iteration reached! Cannot find wα.")
            return
        end
    end
    return wa
end


"Brutal computation of integrand"
function integrand(w,t,α,c,p)

    q = qα(w,t,α,p)
    dqdt = dqαdt(w,t,α,p)
    
    E = 1 + q^2*p[:D] + w^2*p[:F]
    O = -2*q*w*(p[:Σ]^2-p[:N]^2)*cos(p[:ϕ])*sin(p[:ϕ])

    return real(
        ((E^2 + O^2)*Mα(w,q,α,j,p) -
         2*E*O*Nα(w,q,α,j,p))/((E^2-O^2)^2) * dqdt)

end

"regularised integrand, removing singularity at w=Tα"
function integrand_smooth(Tα,χ,t,α,c,p)
    # transformation is the following:
    w = Tα*sin(χ)

    # and now we integrate in χ from 0 to π/2.

    # the derivative term dqdt *dw can thus be replaced by:
    dqdtdw = t*cos(p[:θ]) + 1.0im*sin(p[:θ])*cos(χ)*Tα
    
    q = qα(w,t,α,p)
    
    E = 1 + q^2*p[:D] + w^2*p[:F]
    O = -2*q*w*(p[:Σ]^2-p[:N]^2)*cos(p[:ϕ])*sin(p[:ϕ])

    return real(
        ((E^2 + O^2)*Mα(w,q,α,c,p) -
         2*E*O*Nα(w,q,α,c,p))/((E^2-O^2)^2) * dqdtdw)
end

"fully regularised integrand, with all singularities removed, including that at w=wα"
function integrand_smooth_nonsingular(Tα, χ, wα, t, α, c ,p)
    w = Tα*sin(χ)
    dqdtdw = t*cos(p[:θ]) + 1.0im*sin(p[:θ])*cos(χ)*Tα
    q = qα(w,t,α,p)
    E = 1 + q^2*p[:D] + w^2*p[:F]
    O = -2*q*w*(p[:Σ]^2-p[:N]^2)*cos(p[:ϕ])*sin(p[:ϕ])

    return real(
        ((E^2 + O^2)*Mα(w,q,α,c,p) -
         2*E*O*Nα(w,q,α,c,p))/((E^2-O^2)^2) * dqdtdw) - Sα(w,wα,t,α,c,p)*Tα*cos(χ)
    
end

function XαYα(wα,t,α,c,p)
    q = qα(wα,t,α,p)
    τ = τα(α,p)

    a3p = dqαdw(wα,t,α,p)
    da3p = d2qαdw2(wα,t,α,p)

    a_14 = a14(wα,q,p,a3p)
    a_15 = a15(wα,q,p,a3p,da3p)
    a_18 = t/√(t^2 - wα^2-τ^2)*cos(p[:θ]) + 1.0im*sin(p[:θ])
    a_19 = wα*t*(t^2 - wα^2-τ^2)^(-3/2)*cos(p[:θ])
    
    V = Vα(wα,q,α,c,p)
    W = Wα(wα,q,α,c,p,a3p)
    
    X = a_18*V/(8*a_14^2)
    Y = ((a_19 - a_15*a_18)*V + a_18*W)/(8*a_14^2)

    return (X,Y)
end


function Sα(w,wα,t,α,c,p)
    (X,Y) = XαYα(wα,t,α,c,p)

    return real(
        X/(w-wα)^2 + Y/(w-wα)
    )
end

function integralofSα(Tα,wα,t,α,c,p)
    (X,Y) = XαYα(wα,t,α,c,p)

    return real(
        X*Tα/(wα*(wα-Tα)) + Y*(log(abs(1-Tα/wα)) + 1.0im*angle(1-Tα/wα))
    )
end


function Mα(w,q,α,c,p)
    if α=='d'
        if c=="ux"
            return -q^2*(cos(p[:ϕ]))^2 - w^2*(sin(p[:ϕ]))^2
        elseif c=="uy"
            return (-q^2+w^2)*cos(p[:ϕ])*sin(p[:ϕ])
        elseif c=="uz"
            return -1.0im*q*md(w,q)*cos(p[:ϕ])
        elseif c=="zx"
            return md(w,q)*(q^2*cos(p[:ϕ])^2 + w^2*sin(p[:ϕ])^2)
        elseif c=="zy"
            return md(w,q)*(q^2-w^2)*cos(p[:ϕ])*sin(p[:ϕ])
        elseif c=="zz"
            return 1.0im*q*(q^2 + w^2 + 0.5*p[:cdcs]^2)*cos(p[:ϕ])
        end
    elseif α=='s'
        if c=="ux"
            return q^2*(cos(p[:ϕ]))^2 + w^2*(sin(p[:ϕ]))^2 + 0.5*p[:cdcs]^2
        elseif c=="uy"
            return (q^2-w^2)*cos(p[:ϕ])*sin(p[:ϕ])
        elseif c=="uz"
            m = ms(w,q,p)
            return 1.0im*q*(m - 0.5*p[:cdcs]^2/m)*cos(p[:ϕ])
        elseif c=="zx"
            m = ms(w,q,p)
            return (-q^2*cos(p[:ϕ])^2 - w^2*sin(p[:ϕ])^2)*(m - 0.25*p[:cdcs]^2/m) - 0.25*p[:cdcs]^2*m
        elseif c=="zy"
            m = ms(w,q,p)
            return (-q^2+w^2)*cos(p[:ϕ])*sin(p[:ϕ])*(m - 0.25*p[:cdcs]^2/m)
        elseif c=="zz"
            return -1.0im*q*(q^2 + w^2 + 0.5*p[:cdcs]^2)*cos(p[:ϕ])
        end
    end
end

function Nα(w,q,α,c,p)
    if α=='d'
        if c=="ux"
            return 2*q*w*cos(p[:ϕ])*sin(p[:ϕ])
        elseif c=="uy"
            return -q*w*(cos(p[:ϕ])^2 - sin(p[:ϕ])^2)
        elseif c=="uz"
            return 1.0im*w*md(w,q)*sin(p[:ϕ])
        elseif c=="zx"
            return -2*q*w*md(w,q)*cos(p[:ϕ])*sin(p[:ϕ])
        elseif c=="zy"
            return q*w*md(w,q)*(cos(p[:ϕ])^2 - sin(p[:ϕ])^2)
        elseif c=="zz"
            return -1.0im*w*(q^2 + w^2 + 0.5*p[:cdcs]^2)*sin(p[:ϕ])
        end
    elseif α=='s'
        if c=="ux"
            return -2*q*w*cos(p[:ϕ])*sin(p[:ϕ])
        elseif c=="uy"
            return q*w*(cos(p[:ϕ])^2 - sin(p[:ϕ])^2)
        elseif c=="uz"
            m = ms(w,q,p)
            return -1.0im*w*(m-0.5*p[:cdcs]^2/m)*sin(p[:ϕ])
        elseif c=="zx"
            m = ms(w,q,p)
            return 2*q*w*(m-0.25*p[:cdcs]^2/m)*cos(p[:ϕ])*sin(p[:ϕ])
        elseif c=="zy"
            m = ms(w,q,p)
            return -q*w*(m-0.25*p[:cdcs]^2/m)*(cos(p[:ϕ])^2 - sin(p[:ϕ])^2)
        elseif c=="zz"
            return 1.0im*w*(q^2 + w^2 + 0.5*p[:cdcs]^2)*sin(p[:ϕ])
        end
    end
end

function residue(α,c,w,p)

    q = qσν(w,p)
    
    a_1  = a1(w,p)
    a_2  = a2(w,q,α,p)
    a_11 = a11(w,q,α,p)
    a_12 = a12(w,q,α,p)
    U = Uα(w,q,α,c,p)
    V = Vα(w,q,α,c,p)
    W = Wα(w,q,α,c,p,dqσνdw(w,p))

    return 2π*real(
        1.0im*( -(U + 1.0im*V*p[:D]/a_1)*a_2/(8*a_1^2) +
                a_11*a_2*W + a_12*V)
    )
end

function Uα(w,q,α,c,p)
    if α=='d'
        if c=="ux"
            return -2*a4(w,q,p)*cos(p[:ϕ])
        elseif c=="uy"
            return -a6(w,q,p)
        elseif c=="uz"
            m = md(w,q)
            return -1.0im*(m*cos(p[:ϕ]) + q*a4(w,q,p)/m)
        elseif c=="zx"
            m = md(w,q)
            a_4 = a4(w,q,p)
            return 2*a_4*m*cos(p[:ϕ]) + q*a_4^2/m
        elseif c=="zy"
            m = md(w,q)
            return a6(w,q,p)*m + q*a5(w,q,p)/m
        elseif c=="zz"
            return 1.0im*(a13(w,q,p)*cos(p[:ϕ]) + 2*q*a4(w,q,p))
        end
    elseif α=='s'
        if c=="ux"
            return 2*a4(w,q,p)*cos(p[:ϕ])
        elseif c=="uy"
            return a6(w,q,p)
        elseif c=="uz"
            return 1.0im*(a71(w,q,p)*cos(p[:ϕ]) + a4(w,q,p)*a73(w,q,p))
        elseif c=="zx"
            a_4 = a4(w,q,p)
            return -2*a_4*a72(w,q,p)*cos(p[:ϕ]) - a_4^2*a74(w,q,p) -
                0.25*p[:cdcs]^2*q/a7(w,q,p)
        elseif c=="zy"
            return -a6(w,q,p)*a72(w,q,p) -a5(w,q,p)*a74(w,q,p)
        elseif c=="zz"
            return -1.0im*(a13(w,q,p)*cos(p[:ϕ]) + 2*q*a4(w,q,p))
        end
    end
end

function Vα(w,q,α,c,p)
    if α=='d'
        if c=="ux"
            return -a4(w,q,p)^2
        elseif c=="uy"
            return -a5(w,q,p)
        elseif c=="uz"
            return -1.0im*a4(w,q,p)*√(q^2+w^2+1)
        elseif c=="zx"
            return a4(w,q,p)^2*√(q^2+w^2+1)
        elseif c=="zy"
            return a5(w,q,p)*√(q^2+w^2+1)
        elseif c=="zz"
            return 1.0im*a4(w,q,p)*a13(w,q,p)
        end
    elseif α=='s'
        if c=="ux"
            return a4(w,q,p)^2 + 0.5*p[:cdcs]^2
        elseif c=="uy"
            return a5(w,q,p)
        elseif c=="uz"
            return 1.0im*a4(w,q,p)*a71(w,q,p)
        elseif c=="zx"
            return -a4(w,q,p)^2*a72(w,q,p) - 0.25*p[:cdcs]^2*a7(w,q,p)
        elseif c=="zy"
            return -a5(w,q,p)*a72(w,q,p)
        elseif c=="zz"
            return -1.0im*a4(w,q,p)*a13(w,q,p)
        end
    end
end

function Wα(w,q,α,c,p,a3)
    if α=='d'
        if c=="ux"
            return -2*a4(w,q,p)*a8(a3,p)
        elseif c=="uy"
            return -a9(a3,w,q,p)
        elseif c=="uz"
            m = md(w,q)
            return -1.0im*(a8(a3,p)*m + a4(w,q,p)*a10(a3,w,q)/m)
        elseif c=="zx"
            m = md(w,q)
            return 2*a4(w,q,p)*a8(a3,p)*m +
                a4(w,q,p)^2*a10(a3,w,q)/m
        elseif c=="zy"
            m = md(w,q)
            return a9(a3,w,q,p)*m + a5(w,q,p)*a10(a3,w,q)/m
        elseif c=="zz"
            return 1.0im*(a8(a3,p)*a13(w,q,p) +
                          2*a4(w,q,p)*a10(a3,w,q))
        end
    elseif α=='s'
        if c=="ux"
            return 2*a4(w,q,p)*a8(a3,p)
        elseif c=="uy"
            return a9(a3,w,q,p)
        elseif c=="uz"
            return 1.0im*(a8(a3,p)*a71(w,q,p) +
                          a4(w,q,p)*a10(a3,w,q)*a73(w,q,p)/q)
        elseif c=="zx"
            return -2*a4(w,q,p)*a8(a3,p)*a72(w,q,p)-
            (a4(w,q,p)^2*a74(w,q,p)/q +
             0.25*p[:cdcs]^2/a7(w,q,p))*a10(a3,w,q)
        elseif c=="zy"
            return -a9(a3,w,q,p)*a72(w,q,p) -
                a5(w,q,p)*a10(a3,w,q)*a74(w,q,p)/q
        elseif c=="zz"
            return -1.0im*(a8(a3,p)*a13(w,q,p) +
                           2*a4(w,q,p)*a10(a3,w,q))
        end
    end
end

function ms(w,q,p)
    return √(q^2+w^2+p[:cdcs]^2)
end

function md(w,q)
    return √(q^2+w^2+1)
end

function mα(w,q,α,p)
    if α=='d'
        return md(w,q)
    elseif α=='s'
        return ms(w,q,p)
    end
end

function a1(w,p)
    return √(w^2*p[:Σ]^2*p[:N]^2 + p[:D])
end

function a2(w,q,α,p)
    dqdw = dqσνdw(w,p)
    return 1/(-1.0im*sin(p[:θ])*dqdw +
              cos(p[:θ])*(q*dqdw + w)/mα(w,q,α,p))
end

function ineq(p)
    return ((p[:Σ]^2-p[:N]^2)*sin(2*p[:ϕ])>0)
end

function a4(w,q,p)
    if ineq(p) 
        return q*cos(p[:ϕ]) - w*sin(p[:ϕ])
    else
        return q*cos(p[:ϕ]) + w*sin(p[:ϕ])
    end
end

function a5(w,q,p)
    if ineq(p)
        return (q^2-w^2)*cos(p[:ϕ])*sin(p[:ϕ]) +
            q*w*(cos(p[:ϕ])^2 - sin(p[:ϕ])^2)
    else
        return (q^2-w^2)*cos(p[:ϕ])*sin(p[:ϕ]) -
            q*w*(cos(p[:ϕ])^2 - sin(p[:ϕ])^2)
    end
end

function a6(w,q,p)
    if ineq(p)
        return 2*q*cos(p[:ϕ])*sin(p[:ϕ]) +
            w*(cos(p[:ϕ])^2 - sin(p[:ϕ])^2)
    else
        return 2*q*cos(p[:ϕ])*sin(p[:ϕ]) -
            w*(cos(p[:ϕ])^2 - sin(p[:ϕ])^2)
    end
end

function a7(w,q,p)
    return ms(w,q,p)
end

function a71(w,q,p)
    m = ms(w,q,p)
    return m - 0.5*p[:cdcs]^2/m
end

function a72(w,q,p)
    m = ms(w,q,p)
    return m - 0.25*p[:cdcs]^2/m
end

function a73(w,q,p)
    m = ms(w,q,p)
    return q/m + 0.5*q*p[:cdcs]^2/m^3
end

function a74(w,q,p)
    m = ms(w,q,p)
    return q/m + 0.25*q*p[:cdcs]^2/m^3
end

function a8(a3,p)
    if ineq(p)
        return a3*cos(p[:ϕ]) - sin(p[:ϕ])
    else
        return a3*cos(p[:ϕ]) + sin(p[:ϕ])
    end
end

function a9(a3,w,q,p)
    if ineq(p)
        return 2*(q*a3-w)*cos(p[:ϕ])*sin(p[:ϕ]) +
            (w*a3+q)*(cos(p[:ϕ])^2-sin(p[:ϕ])^2)
    else
        return 2*(q*a3-w)*cos(p[:ϕ])*sin(p[:ϕ]) -
            (w*a3+q)*(cos(p[:ϕ])^2-sin(p[:ϕ])^2)
    end
end

function a10(a3,w,q)
    return q*a3+w
end

function a11(w,q,α,p)
    a3 = dqσνdw(w,p)
    m = mα(w,q,α,p)
    return 1/(8*a1(w,p)^2*(a3 + w*cos(p[:θ])/(q*cos(p[:θ]) - 1.0im*m*sin(p[:θ]))))
end

function a12(w,q,α,p)
    a_11 = a11(w,q,α,p)
    a_1 = a1(w,p)
    a_3 = dqσνdw(w,p)
    m = mα(w,q,α,p)
    da3dw = 1.0im*p[:Σ]^2*p[:N]^2*(1 - w^2*p[:Σ]^2*p[:N]^2/(a_1^2))/(p[:D]*a_1)
    bi = 1/(q*cos(p[:θ]) - 1.0im*m*sin(p[:θ]))
    dmdw = (q*a_3 + w)/m

    return -a2(w,q,α,p)*a_11^2*(
        16*w*p[:Σ]^2*p[:N]^2*(a_3 + w*cos(p[:θ])*bi) + 8*a_1^2*
        (da3dw +
         cos(p[:θ])*bi -
         w*cos(p[:θ])*(a_3*cos(p[:θ]) - 1.0im*sin(p[:θ])*dmdw)*bi^2))
end

function a13(w,q,p)
    return q^2 + w^2 + 0.5*p[:cdcs]^2
end

function a14(w,q,p,dqdw)
    return w*p[:F] - q*p[:Δ] + 1.0*im*a1(w,p)*dqdw
end

function a15(w,q,p,dqdw,d2qdw2)
    a_17 = 1/(w*p[:F] - q*p[:Δ])
    return p[:F]*a_17 - dqdw*(p[:Δ] + q*p[:Σ]^2*p[:N]^2*a_17)*a_17 +
        (1.0*im*a1(w,p)*d2qdw2 - dqdw^2*p[:Σ]^2*p[:N]^2*a_17^2)/a14(w,q,p,dqdw)
end
