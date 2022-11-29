function [tau_sh,tau_lin,cg_sh,cg_lin,k]=tau_dispersion_test()
%% shear current data
Frh       = 0.0;
theta     = 1* pi; 
N         = 127;
h         = 4;
g         = 9.8;
Fr_shear  = 0.1;
alpha     = 10;


%% define functions
dk   = 0.001;
omega_fun_sh   =  @(k) fun_omegaf(N,k,h,theta,Fr_shear,alpha);
group_vec_sh   =  @(k) (omega_fun_sh(k+dk)-omega_fun_sh(k-dk))./(2*dk);

% linear shear profile
% omega_f_lin    =  @(k) sqrt(g*k.*tanh(k*h)+(S1.*tanh(k*h).*cos(theta)/2).^2)-...
%                        S1.*tanh(k*h).*cos(theta)/2;
% group_vec_lin  =  @(k) (omega_f_lin(k+dk)-omega_f_lin(k-dk))./(2*dk);

%% evaluating values
% kh           = 0.01:0.01:20;
k            =  1;  
[omega,cg_Q] = cg_quinn(N,k,h,theta,Fr_shear,alpha);
tau_sh       = omega_fun_sh(k)*sqrt(h/g);
cg_sh        = group_vec_sh(k)/sqrt(g*h)+Fr_shear*cos(theta);

% tau_lin   = omega_f_lin(kh)*sqrt(h/g);
% cg_lin    = group_vec_lin(kh)/sqrt(g*h);

cg_Q      = cg_Q/sqrt(g*h);
tc        = omega./k/sqrt(g*h)+Fr_shear/sqrt(h)*cos(theta);
tc1       = omega_fun_sh(k)/sqrt(g*h)./k+Fr_shear/sqrt(h)*cos(theta);
tc_nvs    = sqrt(tanh(k.*h)./(k.*h))+Fr_shear/sqrt(h)*cos(theta); 
cg_nvs    = sqrt(tanh(k.*h)./(k.*h))/2.*(1+2*(k.*h)./sinh(2*(k.*h)))+Fr_shear/sqrt(h)*cos(theta);

% plot(kh,abs(cg_sh),'--r',kh,abs(cg_Q),'-.b',kh,cg_nvs,'-k')
% plot(kh,abs(abs(cg_sh)-abs(cg_Q))./sqrt(g*h))
plot(k,abs(abs(cg_Q))./abs(cg_nvs),'-.r',k,tc1./tc_nvs,'--r')
% set(gca,'ylim',[0.5 1.5])
end