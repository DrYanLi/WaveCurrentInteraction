function [tc_kc,tc_EL] = obtain_kc_EL_simpsons(I,Ik,z)

%%  inputs
g       = I.g;
h       = I.h;
kv      = Ik.k;
k_z     = Ik.k_i;
kx_z    = Ik.kx_i;
ky_z    = Ik.ky_i;
z_k     = z.z_k;
dz      = z.dz; %h/I.N
   
%%       
th      = tanh(kv.*h);  
c0      = sqrt(g.*th./kv);
   
%% E&L

%% 
sh2sh    = (exp(2*k_z.*z_k)-exp(-4*(1+0.5*z_k/h).*k_z*h))./(1-exp(-4.*k_z*h));
itg_kc   = (kx_z.*I.dUkx_i+ky_z.*I.dUky_i)./k_z.*sh2sh;
itg_kc(isnan(itg_kc)) = 0;

dltL     = simpsons(itg_kc,dz)./c0;
tc_EL    = c0.*(sqrt(1+dltL.^2)-dltL);
tc_kc    = c0-simpsons(itg_kc,dz);
end

