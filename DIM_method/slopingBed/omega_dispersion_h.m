function [tau_sh,tau_lin,cg_sh,cg_lin,k_v]=omega_dispersion_h()
%% shear current data
Frh       = 0.0;
theta     = 1* pi; 
N         = 127;
g         = 9.8;
Fr_shear  = 0.10;
alpha     = 10;

%
k_v       =  0.2:0.02:3; 
hv         = 4:-0.01:1;
[k_i,h_i]      = meshgrid(k_v,hv);

% kh           = 0.01:0.01:20;
omega_v      = 0.*h_i;
cg_Q         = 0.*h_i;
tc_v         = 0.*h_i;


for i = 1:length(hv)
    [omega_v(i,:),cg_Q(i,:)] = cg_quinn(N,k_v,hv(i),theta,Fr_shear,alpha);
end
omega_nvs = sqrt(tanh(k_i.*h_i)*g.*k_i)+k_i.*Fr_shear.*sqrt(g)*cos(theta);
cg_nvs    = sqrt(tanh(k_i.*h_i).*g./(k_i))/2.*(1+2*(k_i.*h_i)./sinh(2*(k_i.*h_i)))+Fr_shear.*sqrt(g)*cos(theta);
% contourf(k_i,h_i,omega_v)
subplot(1,3,1)
contourf(k_i,h_i,omega_v)
subplot(1,3,2)
contourf(k_i,h_i,omega_nvs)
subplot(1,3,3)

mincg =  min(min(omega_v./omega_nvs));
maxcg = max(max(omega_v./omega_nvs));
stepcg = (maxcg-mincg)/9;
hv = mincg:stepcg:maxcg;
[c,h]=contourf(k_i,h_i,omega_v./omega_nvs,hv)
clabel(c,h,'FontName','Times New Roman','FontSize',9)
a=0;
% tc_v      = omega_v/k_v;
% plot(hv,omega_v,'--r',hv,omega_nvs,'-.k')
% tau_sh       = omega_fun_sh(k_v)*sqrt(h/g);
% cg_sh        = group_vec_sh(k_v)/sqrt(g*h)+Fr_shear*cos(theta);

% tau_lin   = omega_f_lin(kh)*sqrt(h/g);
% cg_lin    = group_vec_lin(kh)/sqrt(g*h);

% cg_Q      = cg_Q/sqrt(g*h);
% tc        = omega./k_v;
% % tc1       = omega_fun_sh(k)/sqrt(g*h)./k+Fr_shear/sqrt(h)*cos(theta);
% % tc_nvs    = sqrt(tanh(k.*h)./(k.*h))+Fr_shear/sqrt(h)*cos(theta); 
% % cg_nvs    = sqrt(tanh(k.*h)./(k.*h))/2.*(1+2*(k.*h)./sinh(2*(k.*h)))+Fr_shear/sqrt(h)*cos(theta);
% 
% % plot(kh,abs(cg_sh),'--r',kh,abs(cg_Q),'-.b',kh,cg_nvs,'-k')
% % plot(kh,abs(abs(cg_sh)-abs(cg_Q))./sqrt(g*h))
% plot(k_v,abs(abs(cg_Q))./abs(cg_nvs),'-.r',k_v,tc1./tc_nvs,'--r')
% set(gca,'ylim',[0.5 1.5])
end