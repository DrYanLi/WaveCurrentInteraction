clear all
N    = 15;
theta = 0;

%%
sgn = 1;
[tc_str, R_str, U_str] = obtain_data(N,sgn,theta);
tc_PLA       = tc_str.PLA ;
tc_kc        = tc_str.kc;
tc_EL        = tc_str.EL;
kv           = tc_str.kv;
tc_dimM0_1   = tc_str.dim01;
tc_dimM0_2   = tc_str.dim02;
tc_dimM1_1   = tc_str.dim11;

R01    = R_str.R01;
R02    = R_str.R02;
R11    = R_str.R11;
R01_e  = R_str.R01_e;
R02_e  = R_str.R02_e;
R11_e  = R_str.R11_e;

% shear profile
subplot(3,3,1)
plot(U_str.Uz,U_str.z,'-k')
xlabel('$U(z)/\sqrt{gh}$','Interpreter','Latex')
ylabel('$z/h$','Interpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',10)

% tc profile
subplot(3,3,4)
semilogx(kv,tc_kc./tc_PLA,'-.b',kv,tc_EL./tc_PLA,'--r',kv,tc_dimM0_1./tc_PLA,'-..k',kv,tc_dimM0_2./tc_PLA,'-k',kv,tc_dimM1_1./tc_PLA,':k')
legend('K&C_{1st}','E&L_{1st}','M_0 = 1','M_0 = 2','M_{kc0}=1')
set(gca,'xlim',[0.01 100],'FontName','Times New Roman','FontSize',10,'xtick',[0.01 0.1 1 10 100])
xlabel('$kh$','Interpreter','Latex')
ylabel('$\tilde{c}/\tilde{c}_{N}$','Interpreter','latex')

% error estimates
subplot(3,3,7)
semilogx(kv,R01,'-k',kv,R02,'-.k',kv,R11,':k')
legend('M_0 = 1','M_0 = 2','M_{kc0}=1')
set(gca,'xlim',[0.01 100],'FontName','Times New Roman','FontSize',10,'xtick',[0.01 0.1 1 10 100])
xlabel('$kh$','Interpreter','Latex')
ylabel('$\tilde{c}/\tilde{c}_{N}$','Interpreter','latex')
%%

%%
sgn = 2;
[tc_str, R_str, U_str] = obtain_data(N,sgn,theta);
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

% shear profile
subplot(3,3,2)
plot(U_str.Uz,U_str.z,'-k')
xlabel('$U(z)/\sqrt{gh}$','Interpreter','Latex')
ylabel('$z/h$','Interpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',10)

% tc profile
subplot(3,3,5)
semilogx(kv,tc_kc./tc_PLA,'-.b',kv,tc_EL./tc_PLA,'--r',kv,tc_dimM0_1./tc_PLA,'-..k',kv,tc_dimM0_2./tc_PLA,'-k',kv,tc_dimM1_1./tc_PLA,':k')
legend('K&C_{1st}','E&L_{1st}','M_0 = 1','M_0 = 2','M_{kc0}=1')
set(gca,'xlim',[0.01 100],'FontName','Times New Roman','FontSize',10,'xtick',[0.01 0.1 1 10 100])
xlabel('$kh$','Interpreter','Latex')
ylabel('$\tilde{c}/\tilde{c}_{N}$','Interpreter','latex')

% error estimates
subplot(3,3,8)
semilogx(kv,R01,'-k',kv,R02,'-.k',kv,R11,':k')
legend('M_0 = 1','M_0 = 2','M_{kc0}=1')
set(gca,'xlim',[0.01 100],'FontName','Times New Roman','FontSize',10,'xtick',[0.01 0.1 1 10 100])
xlabel('$kh$','Interpreter','Latex')
ylabel('$\tilde{c}/\tilde{c}_{N}$','Interpreter','latex')

%%
sgn = 3;
[tc_str, R_str, U_str] = obtain_data(N,sgn,theta);
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

% shear profile
subplot(3,3,3)
plot(U_str.Uz,U_str.z,'-k')
xlabel('$U(z)/\sqrt{gh}$','Interpreter','Latex')
ylabel('$z/h$','Interpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',10)

% tc profile
subplot(3,3,6)
semilogx(kv,tc_kc./tc_PLA,'-.b',kv,tc_EL./tc_PLA,'--r',kv,tc_dimM0_1./tc_PLA,'-..k',kv,tc_dimM0_2./tc_PLA,'-k',kv,tc_dimM1_1./tc_PLA,':k')
legend('K&C_{1st}','E&L_{1st}','M_0 = 1','M_0 = 2','M_{kc0}=1')
set(gca,'xlim',[0.01 100],'FontName','Times New Roman','FontSize',10,'xtick',[0.01 0.1 1 10 100])
xlabel('$kh$','Interpreter','Latex')
ylabel('$\tilde{c}/\tilde{c}_{N}$','Interpreter','latex')

% error estimates
subplot(3,3,9)
semilogx(kv,R01,'-k',kv,R02,'-.k',kv,R11,':k')
legend('M_0 = 1','M_0 = 2','M_{kc0}=1')
set(gca,'xlim',[0.01 100],'FontName','Times New Roman','FontSize',10,'xtick',[0.01 0.1 1 10 100])
xlabel('$kh$','Interpreter','Latex')
ylabel('$\tilde{c}/\tilde{c}_{N}$','Interpreter','latex')