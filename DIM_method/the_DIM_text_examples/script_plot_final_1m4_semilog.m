clear all
N     = 15; % vertical mesh number
theta = 0;

%%
sgn = 1; % 
[tc_str, R_str, U_str1] = obtain_data(N,sgn,theta);
% tc_PLA       = tc_str.PLA ;
% tc_kc        = tc_str.kc;
% tc_EL        = tc_str.EL;
% kv           = tc_str.kv;
% tc_dimM0_1   = tc_str.dim01;
% tc_dimM0_2   = tc_str.dim02;
% tc_dimM1_1   = tc_str.dim11;
kv           = tc_str.kv;

R01    = R_str.R01;
R02    = R_str.R02;
R11    = R_str.R11;
R01_e  = R_str.R01_e;
R02_e  = R_str.R02_e;
R11_e  = R_str.R11_e;
R_kc   = R_str.R_kc;
R_EL   = R_str.R_EL;

% shear profile

% annotation('textbox',...
%     [0  0.3 0 0],...
%     'String',{'(a)'},...
%     'FontSize',10,...
%     'Edgecolor','none')

% error estimates
subplot(1,4,2)
semilogx(kv,R_kc,'--r',kv,R_EL,'-.b',kv,R01,'-k',kv,R02,'-..k',kv,R11,':k')
legend('K&C_{1st}','E&L_{1st}','M_0 = 1','M_0 = 2','M_{kc0}=1')
set(gca,'xlim',[0.01 100],'FontName','Times New Roman','FontSize',10,'xtick',[0.01 0.1 1 10 100])
xlabel('$kh$','Interpreter','Latex')
ylabel('$R$','Interpreter','latex')
title('(b)','FontName','Times New Roman','FontSize',10)


%%
sgn = 2;
[tc_str, R_str, U_str2] = obtain_data(N,sgn,theta);
tc_PLA       = tc_str.PLA ;
tc_kc        = tc_str.kc;
tc_EL        = tc_str.EL;
tc_kv        = tc_str.kv;
tc_dimM0_1   = tc_str.dim01;
tc_dimM0_2   = tc_str.dim02;
tc_dimM1_1   = tc_str.dim11;

R01    = R_str.R01;
R02    = R_str.R02;
R11    = R_str.R11;
R01_e  = R_str.R01_e;
R02_e  = R_str.R02_e;
R11_e  = R_str.R11_e;
R_kc   = R_str.R_kc;
R_EL   = R_str.R_EL;

% error estimates
subplot(1,4,3)
semilogx(kv,R_kc,'--r',kv,R_EL,'-.b',kv,R01,'-k',kv,R02,'-..k',kv,R11,':k')
legend('K&C_{1st}','E&L_{1st}','M_0 = 1','M_0 = 2','M_{kc0}=1')
set(gca,'xlim',[0.01 100],'FontName','Times New Roman','FontSize',10,'xtick',[0.01 0.1 1 10 100])
xlabel('$kh$','Interpreter','Latex')
ylabel('$R$','Interpreter','latex')
title('(c)','FontName','Times New Roman','FontSize',10)

%%
sgn = 3;
[tc_str, R_str, U_str3] = obtain_data(N,sgn,theta);
tc_PLA       = tc_str.PLA ;
tc_kc        = tc_str.kc;
tc_EL        = tc_str.EL;
tc_kv        = tc_str.kv;
tc_dimM0_1   = tc_str.dim01;
tc_dimM0_2   = tc_str.dim02;
tc_dimM1_1   = tc_str.dim11;

R01    = R_str.R01;
R02    = R_str.R02;
R11    = R_str.R11;
R01_e  = R_str.R01_e;
R02_e  = R_str.R02_e;
R11_e  = R_str.R11_e;
R_kc   = R_str.R_kc;
R_EL   = R_str.R_EL;

% error estimates
subplot(1,4,4)
semilogx(kv,R_kc,'--r',kv,R_EL,'-.b',kv,R01,'-k',kv,R02,'-..k',kv,R11,':k')
legend('K&C_{1st}','E&L_{1st}','M_0 = 1','M_0 = 2','M_{kc0}=1')
set(gca,'xlim',[0.01 100],'FontName','Times New Roman','FontSize',10,'xtick',[0.01 0.1 1 10 100])
xlabel('$kh$','Interpreter','Latex')
ylabel('$R$','Interpreter','latex')
title('(d)','FontName','Times New Roman','FontSize',10)

%%
subplot(1,4,1)
plot(U_str1.Uz,U_str1.z,'-k',U_str2.Uz,U_str2.z,'-.k',U_str3.Uz,U_str3.z,'--k')
legend('Profile 1','Profile 2','Profile 3')
xlabel('$U(z)/\sqrt{gh}$','Interpreter','Latex')
ylabel('$z/h$','Interpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',10,'xlim',[-0.2 0.5],'xtick',[-0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5])
title('(a) shear profiles','FontName','Times New Roman','FontSize',10)
