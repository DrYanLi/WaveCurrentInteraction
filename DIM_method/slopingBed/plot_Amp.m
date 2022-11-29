function [a]=plot_Amp()
%% shear current data
Ins.Frh       = 0.0;

Ins.N         = 127;
Ins.g         = 9.8;
Ins.Fr_shear  = 0.10;
Ins.alpha     = 6;
h_1           = 10;
%

k_v           = 10.^(-1.5:0.01:0.8);
R_A           = 0.*k_v;
RA_nvs        = 0.*k_v;
omega_v       = 0.*k_v;
cg_Q1         = 0.*k_v;   
k_2           = 0.*k_v;
cg_nvs_1      = 0.*k_v;
k_2nvs        = 0.*k_v;
k_1nvs        = 0.*k_v;
%% H = 10
R_A_pi        = 0.*k_v;
for i=1:length(k_v)
    Ins.theta     = 0*pi; 
    [cg_Q1(i),cg_Q2,k_2(i),cg_nvs_1(i),cg_nvs_2,k_1nvs(i),k_2nvs(i),omega_v(i)] = dispersion_ko(k_v(i),h_1,Ins);
    R_A(i)    = sqrt(cg_Q1(i)/cg_Q2);    %A_2/A_1
    RA_nvs(i) = sqrt(cg_nvs_1(i)/cg_nvs_2);
    
    Ins.theta = pi;
    [cg_Q1_x1,cg_Q2_x2,~,~,~,~,~,~] = dispersion_ko(k_v(i),h_1,Ins);
    R_A_pi(i) = sqrt(cg_Q1_x1/cg_Q2_x2);
end
subplot(2,3,1)
plot(omega_v/sqrt(9.8),R_A./RA_nvs,'--r',omega_v/sqrt(9.8),R_A_pi./RA_nvs,'-b')
legend('\theta=0','\theta=\pi')

ylabel('A_2 A_{nvs1}/(A_1A_{nvs2})')
xlabel('$\omega\sqrt{h_2/g}$','Interpreter','Latex')
title('(a) H=10','fontname','times new roman')
set(gca,'fontname','times new roman')

%%
 h_1 = 5;
for i=1:length(k_v)
    Ins.theta     = 0*pi; 
    [cg_Q1(i),cg_Q2,k_2(i),cg_nvs_1(i),cg_nvs_2,k_1nvs,k_2nvs(i),omega_v(i)] = dispersion_ko(k_v(i),h_1,Ins);
    R_A(i)    = sqrt(cg_Q1(i)/cg_Q2);    %A_2/A_1
    RA_nvs(i) = sqrt(cg_nvs_1(i)/cg_nvs_2);
    
    Ins.theta = pi;
    [cg_Q1_x1,cg_Q2_x2,~,~,~,~,~,~] = dispersion_ko(k_v(i),h_1,Ins);
    R_A_pi(i) = sqrt(cg_Q1_x1/cg_Q2_x2);
end
subplot(2,3,2)
% plot(omega_v/sqrt(9.8),R_A,'--r',omega_v/sqrt(9.8),RA_nvs,'-k',omega_v/sqrt(9.8),R_A_pi,'-.r')
plot(omega_v/sqrt(9.8),R_A./RA_nvs,'--r',omega_v/sqrt(9.8),R_A_pi./RA_nvs,'-b')
ylabel('A_2 A_{nvs1}/(A_1A_{nvs2})')
legend('\theta=0','\theta=\pi')
% legend('\theta=0','no vertical shear','\theta=\pi')
% ylabel('A(x_2)/A(x_1)')
xlabel('$\omega\sqrt{h_2/g}$','Interpreter','Latex')
set(gca,'fontname','times new roman')
title('(b) H=5','fontname','times new roman')
%%
h_1 = 2;
for i=1:length(k_v)
    Ins.theta     = 0*pi; 
    [cg_Q1(i),cg_Q2,k_2(i),cg_nvs_1(i),cg_nvs_2,k_1nvs,k_2nvs(i),omega_v(i)] = dispersion_ko(k_v(i),h_1,Ins);
    R_A(i)    = sqrt(cg_Q1(i)/cg_Q2);    %A_2/A_1
    RA_nvs(i) = sqrt(cg_nvs_1(i)/cg_nvs_2);
    
    Ins.theta = pi;
    [cg_Q1_x1,cg_Q2_x2,~,~,~,~,~,~] = dispersion_ko(k_v(i),h_1,Ins);
    R_A_pi(i) = sqrt(cg_Q1_x1/cg_Q2_x2);
end
subplot(2,3,3)
% plot(omega_v/sqrt(9.8),R_A,'--r',omega_v/sqrt(9.8),RA_nvs,'-k',omega_v/sqrt(9.8),R_A_pi,'-.r')
plot(omega_v/sqrt(9.8),R_A./RA_nvs,'--r',omega_v/sqrt(9.8),R_A_pi./RA_nvs,'-b')
ylabel('A_2 A_{nvs1}/(A_1A_{nvs2})')
legend('\theta=0','\theta=\pi')
% legend('\theta=0','no vertical shear','\theta=\pi')
% ylabel('A(x_2)/A(x_1)')
xlabel('$\omega\sqrt{h_2/g}$','Interpreter','Latex')
set(gca,'fontname','times new roman')
title('(c) H=2','fontname','times new roman')
% subplot(3,1,2)
% plot(omega_v,cg_Q1,'--r',omega_v,cg_nvs_1,'-.b')
% 
% subplot(3,1,3)
% plot(omega_v,k_2,'--r',omega_v,k_2nvs,'-.b',omega_v,k_v,'-k')

%%
Ins.Fr_shear  = 0.3;
for i=1:length(k_v)
    Ins.theta     = 0*pi; 
    [cg_Q1(i),cg_Q2,k_2(i),cg_nvs_1(i),cg_nvs_2,k_1nvs(i),k_2nvs(i),omega_v(i)] = dispersion_ko(k_v(i),h_1,Ins);
    R_A(i)    = sqrt(cg_Q1(i)/cg_Q2);    %A_2/A_1
    RA_nvs(i) = sqrt(cg_nvs_1(i)/cg_nvs_2);
    
    Ins.theta = pi;
    [cg_Q1_x1,cg_Q2_x2,~,~,~,~,~,~] = dispersion_ko(k_v(i),h_1,Ins);
    R_A_pi(i) = sqrt(cg_Q1_x1/cg_Q2_x2);
end
subplot(2,3,4)
% plot(omega_v/sqrt(9.8),R_A,'--r',omega_v/sqrt(9.8),RA_nvs,'-k',omega_v/sqrt(9.8),R_A_pi,'-.r')
plot(omega_v/sqrt(9.8),R_A./RA_nvs,'--r',omega_v/sqrt(9.8),R_A_pi./RA_nvs,'-b')
% plot(k_2,R_A./RA_nvs,'--r',k_2,R_A_pi./RA_nvs,'-b')
legend('\theta=0','\theta=\pi')
% legend('\theta=0','no vertical shear','\theta=\pi')
ylabel('A_2 A_{nvs1}/(A_1A_{nvs2})')
xlabel('$\omega\sqrt{h_2/g}$','Interpreter','Latex')
title('(d) H=10','fontname','times new roman')
set(gca,'fontname','times new roman')

%%
 h_1 = 5;
for i=1:length(k_v)
    Ins.theta     = 0*pi; 
    [cg_Q1(i),cg_Q2,k_2(i),cg_nvs_1(i),cg_nvs_2,k_1nvs,k_2nvs(i),omega_v(i)] = dispersion_ko(k_v(i),h_1,Ins);
    R_A(i)    = sqrt(cg_Q1(i)/cg_Q2);    %A_2/A_1
    RA_nvs(i) = sqrt(cg_nvs_1(i)/cg_nvs_2);
    
    Ins.theta = pi;
    [cg_Q1_x1,cg_Q2_x2,~,~,~,~,~,~] = dispersion_ko(k_v(i),h_1,Ins);
    R_A_pi(i) = sqrt(cg_Q1_x1/cg_Q2_x2);
end
subplot(2,3,5)
% plot(omega_v/sqrt(9.8),R_A,'--r',omega_v/sqrt(9.8),RA_nvs,'-k',omega_v/sqrt(9.8),R_A_pi,'-.r')
plot(omega_v/sqrt(9.8),R_A./RA_nvs,'--r',omega_v/sqrt(9.8),R_A_pi./RA_nvs,'-b')
ylabel('A_2 A_{nvs1}/(A_1A_{nvs2})')
legend('\theta=0','\theta=\pi')
% legend('\theta=0','no vertical shear','\theta=\pi')
% ylabel('A(x_2)/A(x_1)')
xlabel('$\omega\sqrt{h_2/g}$','Interpreter','Latex')
set(gca,'fontname','times new roman')
title('(e) H=5','fontname','times new roman')
%%
h_1 = 2;
for i=1:length(k_v)
    Ins.theta     = 0*pi; 
    [cg_Q1(i),cg_Q2,k_2(i),cg_nvs_1(i),cg_nvs_2,k_1nvs,k_2nvs(i),omega_v(i)] = dispersion_ko(k_v(i),h_1,Ins);
    R_A(i)    = sqrt(cg_Q1(i)/cg_Q2);    %A_2/A_1
    RA_nvs(i) = sqrt(cg_nvs_1(i)/cg_nvs_2);
    
    Ins.theta = pi;
    [cg_Q1_x1,cg_Q2_x2,~,~,~,~,~,~] = dispersion_ko(k_v(i),h_1,Ins);
    R_A_pi(i) = sqrt(cg_Q1_x1/cg_Q2_x2);
end
subplot(2,3,6)
% plot(omega_v/sqrt(9.8),R_A,'--r',omega_v/sqrt(9.8),RA_nvs,'-k',omega_v/sqrt(9.8),R_A_pi,'-.r')
plot(omega_v/sqrt(9.8),R_A./RA_nvs,'--r',omega_v/sqrt(9.8),R_A_pi./RA_nvs,'-b')
ylabel('A_2 A_{nvs1}/(A_1A_{nvs2})')
legend('\theta=0','\theta=\pi')
% legend('\theta=0','no vertical shear','\theta=\pi')
% ylabel('A(x_2)/A(x_1)')
xlabel('$\omega\sqrt{h_2/g}$','Interpreter','Latex')
set(gca,'fontname','times new roman')
title('(f) H=2','fontname','times new roman')
end

function [cg_Q1,cg_Q2,k_2,cg_nvs_1,cg_nvs_2,k_1nvs,k_2nvs,omega_v] = dispersion_ko(k_1,h_1,Ins)
% h_2 is set zero by default
Frh       = Ins.Frh;
theta     = Ins.theta;
N         = Ins.N;
g         = Ins.g;
Fr_shear  = Ins.Fr_shear;
alpha     = Ins.alpha;

Fr_0      = 0;% if it equals zero, then we define the velocities in relative to the surface velocity
%
[omega_v,cg_Q1] = cg_quinn(N,k_1,h_1,theta,Fr_shear,alpha);
omega_v         = omega_v + k_1.*Fr_0.*sqrt(g)*cos(theta);
%% with vertical shear
h_2        = 1;
omega_fk   = @(k) omega_fun(N,k,h_2,theta,Fr_shear,alpha)+ k.*Fr_0.*sqrt(g)*cos(theta);
omega_funf = @(k) omega_fk(k)-omega_v; 

k_int            = 1.2*k_1;
% k0  = fzero(group_vel,kint);
k_2              = fzero(omega_funf,k_int);
[omega_v2,cg_Q2] = cg_quinn(N,k_2,h_2,theta,Fr_shear,alpha);

%% no vertial shear

omega_fnvs_1   = @(k) sqrt(tanh(k.*h_1)*g.*k)+k.*Fr_0.*sqrt(g)*cos(theta)-omega_v;
omega_fnvs_2   = @(k) sqrt(tanh(k.*h_2)*g.*k)+k.*Fr_0.*sqrt(g)*cos(theta)-omega_v;

k_1nvs         = fzero(omega_fnvs_1,k_1);
k_2nvs         = fzero(omega_fnvs_2,k_2);


cg_nvs_1       = sqrt(tanh(k_1nvs.*h_1).*g./(k_1nvs))/2.*(1+2*(k_1nvs.*h_1)./sinh(2*(k_1nvs.*h_1)))+Fr_0.*sqrt(g)*cos(theta);
cg_nvs_2       = sqrt(tanh(k_2nvs.*h_2).*g./(k_2nvs))/2.*(1+2*(k_2nvs.*h_2)./sinh(2*(k_2nvs.*h_2)))+Fr_0.*sqrt(g)*cos(theta);
end