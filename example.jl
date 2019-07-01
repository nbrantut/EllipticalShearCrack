## header
using PyPlot, BaseTools

include("ell_shearcrack.jl")

## parameters

# approx rayleight speed
cRcs(ν) = (0.862 + 1.14ν)/(1+ν)

t = range(0, stop=5, length=1500)
θ = deg2rad(85)
ϕ = deg2rad(89)
κ = 1.0
cdcs = √3
ν=.25
(Σ,N) = (κ*cRcs(ν)/cdcs, κ/cdcs)

ux = ellipticalshearcrack(t, "ux", θ, ϕ, Σ, N, cdcs)
uy = ellipticalshearcrack(t, "uy", θ, ϕ, Σ, N, cdcs)
uz = ellipticalshearcrack(t, "uz", θ, ϕ, Σ, N, cdcs)
zx = ellipticalshearcrack(t, "zx", θ, ϕ, Σ, N, cdcs)
zy = ellipticalshearcrack(t, "zy", θ, ϕ, Σ, N, cdcs)
zz = ellipticalshearcrack(t, "zz", θ, ϕ, Σ, N, cdcs)

## plots
subplot(4,1,1)
plot(t,ux,"k")
subplot(4,1,2)
plot(t, uz, "k")
subplot(4,1,3)
plot(t, zx, "k")
subplot(4,1,4)
plot(t, zz, "k")
