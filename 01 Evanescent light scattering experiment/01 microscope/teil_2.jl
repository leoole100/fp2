# programm funktioniert möglicherweise, aber es fehlt noch an brauchbaren test-daten
# aktuell wirft es einen fehler, weil die zu optimierende fehlerfunktion im suchintervall kein minimum hat

using Roots, ForwardDiff

k_B = 1
eta = 1
beta = 1

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
    sigmaq/(2*delta_t)
end 

function procedure(I_z,range_for_I_0,delta_t,R,T)
    D_0 = calc_D_0(R,T)
    t_steps = [n*delta_t for n in range(1,round(Int,(length(I_z)-1)/delta_t)-1)]
    err_func(I_0) = sum((calc_D_z.(sigmaq.(slice_dz.(t_steps,delta_t,Ref(delta_z(find_z.(beta,I_z,I_0))))),delta_t).-D_orth.(D_0,R,sum.(slice_dz.(t_steps,delta_t,Ref(find_z.(beta,I_z,I_0))))./delta_t)).^2)
    D(f)=x -> ForwardDiff.derivative(err_func,float(x))
    find_zero(D(err_func),range_for_I_0)
end


I_z = exp.(randn(50))
range_for_I_0 = (1,20)
delta_t = 10
R = 1
T = 100

procedure(I_z,range_for_I_0,delta_t,R,T)