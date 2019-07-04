using PyPlot

include("ell_shearcrack.jl")

# approx rayleight speed
cRcs(ν) = (0.862 + 1.14ν)/(1+ν)

t = range(0, stop=5, length=1500)
θ = deg2rad(85)
ϕ = deg2rad(89)
cdcs = √3
ν=.25

figure()

for κ in [0.5;0.7;0.9;1.0]
    (Σ,N) = (κ*cRcs(ν)/cdcs, κ/cdcs)
    ux = ellipticalshearcrack(t, "ux", θ, ϕ, Σ, N, cdcs;rtol=1e-5)
    plot(t, ux, linewidth=1.0)
end

plot(1,0,"k.")
plot(cdcs,0,"k.")
annotate("P",xy=[1,5], ha="right")
annotate("S",xy=[cdcs,5], ha="right")
xlim(0,5)
ylim(-75,150)
xlabel("normalised time")
ylabel("normalised x-acceleration")
legend(["κ=0.5";"κ=0.7";"κ=0.9";"κ=1.0"], loc="best")
