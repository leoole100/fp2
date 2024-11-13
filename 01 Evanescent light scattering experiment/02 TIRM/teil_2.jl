# programm funktioniert möglicherweise, aber es fehlt noch an brauchbaren test-daten
# aktuell wirft es einen fehler, weil die zu optimierende fehlerfunktion im suchintervall kein minimum hat

using Roots: fzero
using ForwardDiff: derivative
using CairoMakie
using Optim: optimize, LBFGS
using LsqFit: curve_fit
using Glob: glob
using DataFrames: DataFrame
using KernelDensity: kde
import CSV, Images
import Optim
include("../functions.jl")

cd(@__DIR__)

# %%
k_B = 1.38e-23
eta = 1.0019e-3
beta = 1/10e-9
TIMESTEP = 1e-3

function read_data()
    path = "../data/TIRM/"
    paths = glob(path * "*/*.dat")
    paths = filter(p -> occursin("failed", lowercase(p)) == false, paths)
    Ia = []
    pa = []
    for p in paths
        I = DataFrame(CSV.File(p, header=["t", "I"])).I
        push!(Ia, I)
    push!(pa, split(p, "/")[end-1])
    end
    df = DataFrame(I=Ia, p=pa)
    df.l = map(p -> split(p, " ")[3]*" "*split(p, " ")[2], df[:,:p])
    df.ot = map(p -> parse(Float64, split(p, " ")[3][3:end]), df[:,:p])
    df
end

function calc_D_0(R,T)
    k_B*T/(6*pi*eta*R)
end

function D_orth(D_0, R, z)
    D_0/(R/z+0.2*log(R/z)+0.9712)
end

function slice_dz(t_0,delta_t,dz)
    dz[t_0:min(t_0+delta_t,length(dz))]
end

function find_z(beta, I_z, I_0)
    max((log(I_0) - log(I_z))/beta,0)
end

function delta_z(z_traj)
    (circshift(z_traj,-1)-z_traj)[1:length(z_traj)-1]
end

function sigmaq(delta_z)
    sum((delta_z .- sum(delta_z)/length(delta_z)).^2/(length(delta_z)-1))
end

function calc_D_z(sigmaq,delta_t)
    sigmaq/(2*delta_t*TIMESTEP)
end 

function estimate_I0(I_z,guessI0,delta_t,R,T)
    D_0 = calc_D_0(R,T)
    t_steps = [n*delta_t for n in range(1,round(Int,(length(I_z)-1)/delta_t)-1)]
    err_func(I_0) = sum((calc_D_z.(sigmaq.(slice_dz.(t_steps,delta_t,Ref(delta_z(find_z.(beta,I_z,I_0))))),delta_t).-D_orth.(D_0,R,sum.(slice_dz.(t_steps,delta_t,Ref(find_z.(beta,I_z,I_0))))./delta_t)).^2)
    D(f)=x -> derivative(err_func,float(x))
    I_0 = fzero(D(err_func), guessI0)
    println(err_func(I_0))
    I_0
end
#%%

data = read_data()
#data.I[1]
data.ot

I_z = data.I[2]
#I_z = exp.(beta * randn(Integer(1e6)))
guessI0 = 10
delta_t = 1000
R = 1.86e-6
T = 273+22

estimate_I0(I_z,guessI0, delta_t, R, T)

#%%
D_0 = calc_D_0(R,T)
t_steps = [n*delta_t for n in range(1,round(Int,(length(I_z)-1)/delta_t)-1)]
model(p) = sum((calc_D_z.(sigmaq.(slice_dz.(t_steps,delta_t,Ref(delta_z(find_z.(p[1],I_z,p[2]))))),delta_t).-D_orth.(D_0,R,sum.(slice_dz.(t_steps,delta_t,Ref(find_z.(p[1],I_z,p[2]))))./delta_t)).^2)
#%%
function find_opt(err_func,p0,num_iter)
    params = copy(p0)
    dp = zeros(length(params))
    for iter in range(1,num_iter)
            #i = 1+Int(round(length(params)*rand()-0.5))
            i = ((length(params)+1)*iter)%length(params)+1
            er = err_func(params)
            if(er<1e-2)
                return params
            end
            #println(er)
            dp[i] = 1e-1*params[i]
            eru = err_func(params.+dp)
            erl = err_func(params.-dp)
            if erl!=eru
                dp[i] = dp[i]*(-2e-1*er/(eru-erl))
                if err_func(params.+dp)<=er
                    params .+= dp
                end
                dp[i] = 0
            end
    end
    params
end
(vrange,trahi) = dist(I_z[1:Int(9e5)])
mask = (maximum(trahi).-trahi).<0.1
model(v,p) = p[4]/(sqrt(2*pi)*p[1]).*exp.(-p[3]^2*(log.(abs.(1 .-v)).-log.(abs(p[2]))).^2 ./(2*p[1]^2))
lb = [1e-9,minimum(I_z),10e-9,0]
ub = [1e-7,maximum(I_z),1e-1,10000]
# p0 = [2.35e-8,0.53,1e-7,1.83e-7]
p0 = [2.35e-8,1-sum(vrange[mask])/sum(mask),1e-7,2e-7]
#p0 = [2.35e-8,0.53,1e-7,2e-7] # sigma, I(z_0), beta, Normierungskonstante
#err_func(p) = sum((model(vrange,p).-trahi).^2)
err_func(p) = 1e-3*sum((model(vrange,p).-trahi).*vrange)^2+5e0*(maximum(model(vrange,p))-maximum(trahi))^4+3e-4*(sum(model(vrange,p))-sum(trahi))^2
#pf = optimize(err_func,p0, Optim.BFGS(),Optim.Options(
#    iterations=100, 
#    time_limit=1,
#))
params = find_opt(err_func,p0,1000)#Optim.minimizer(pf)
I_z0 = params[2]#Optim.minimizer(pf)[2]
beta = params[3]#Optim.minimizer(pf)[3]
f = Figure()
i = Axis(f[1,1])
lines!(vrange,trahi,label="trahi")
lines!(vrange,model(vrange,p0),label="p0")
lines!(vrange,model(vrange,params),label="pf",linestyle=:dot)
axislegend()
save("../figures/02_teil2_01_fit.pdf", f)
f

#%%
#copy of part above, to fit parts of trajectory

z = (log.(I_z).-log(I_z0)).*-beta

num_steps=Int(9e4)
sigma_f = zeros(Int(9e5/num_steps))
I_z0_f = zeros(Int(9e5/num_steps))
#beta_f = zeros(Int(9e5/num_steps))
amp = zeros(Int(9e5/num_steps))
z_loc = zeros(Int(9e5/num_steps))

mask = (maximum(trahi).-trahi).<0.1
fak = 1+1e-3*log(num_steps/9e5)
p_start = [fak*2.35e-8,1-sum(vrange[mask])/sum(mask),fak*2e-7]
params = copy(p_start)

for i in range(1,Int(9e5/num_steps))
    z_loc[i] = sum(z[(i-1)*num_steps+1:i*num_steps])/num_steps
    (vrange,trahi) = dist(I_z[(i-1)*num_steps+1:i*num_steps])
    model(v,p) = p[3]/(sqrt(2*pi)*p[1]).*exp.(-beta^2*(log.(abs.(1 .-v)).-log.(abs(p[2]))).^2 ./(2*p[1]^2))
    #lb = [1e-9,minimum(I_z),10e-9,0]
    #ub = [1e-7,maximum(I_z),1e-1,10000]
    # p0 = [2.35e-8,0.53,1e-7,1.83e-7]
    #p0 = [2.35e-8,0.53,1e-7,2e-7] # sigma, I(z_0), beta, Normierungskonstante
    #mask = (maximum(trahi).-trahi).<0.1
    #fak = 1+1e-3*log(num_steps/9e5)
    #p0 = [fak*2.35e-8,1-sum(vrange[mask])/sum(mask),fak*2e-7]
    err_func(p) = 1e-3*sum((model(vrange,p).-trahi).*vrange)^2+5e0*(maximum(model(vrange,p))-maximum(trahi))^2+3e-4*(sum(model(vrange,p))-sum(trahi))^2
    #err_func(p) = sum((model(vrange,p).-trahi).^4)
    #pf = optimize(err_func,p0, Optim.BFGS(),Optim.Options(
    #    iterations=100, 
    #    time_limit=1,
    #))
    p0 = p_start
    if err_func(params)<err_func(p_start)
        p0 = params
    end
    params = find_opt(err_func,p0,100)#Optim.minimizer(pf)
    (sigma_f[i],I_z0_f[i],amp[i]) = params
    println("Step ", i, ", σ=", params[1], " I(z_0)=", params[2], " Error: ", err_func(params))
end
D_f = sigma_f.^2/(2*num_steps*TIMESTEP)
D_0 = calc_D_0(R,T)

f = Figure()
i = Axis(f[1,1])
i2 = Axis(f[2,1])
scatter!(z_loc,D_f,label="D über z")
#lines!(z_loc,D_orth.(D_0,R,z_loc),linestyle=:dot)
lines!(i,sigma_f,label="σ")
#lines!(i2,z,label="z")
#lines!(vrange,trahi,label="trahi")
#lines!(vrange,model(vrange,p0),label="p0")
#lines!(vrange,model(vrange,Optim.minimizer(pf)),label="pf",linestyle=:dot)
axislegend()
save("../figures/02_teil2_01_fit.pdf", f)
f


#%%
f = Figure()
i = Axis(f[1,1], ylabel="intensity")
ih = Axis(f[1,2])
lines!(i, I_z)
density!(ih, I_z, direction=:y)
f
#%%

# I0 = estimate_I0(I_z,guessI0,delta_t,R,T)
I0 = maximum(I_z)


f = Figure()
i = Axis(f[1,1], ylabel="intensity")
ih = Axis(f[1,2])
z = Axis(f[2,1], xlabel="time", ylabel="z")
zh = Axis(f[2,2])
linkxaxes!(i, z)
linkyaxes!(z, zh)
linkyaxes!(i, ih)
hidexdecorations!(i, grid=false)
hidedecorations!(ih, grid=false)
hidedecorations!(zh, grid=false)
colsize!(f.layout, 1, Auto(3))
xlims!(zh, 0, nothing)
lines!(i, I_z)
density!(ih, I_z, direction=:y)
lines!(z, find_z.(beta,I_z,I0))
density!(zh, find_z.(beta,I_z,I0), direction=:y)
save("../figures/02_path.png", f)
f