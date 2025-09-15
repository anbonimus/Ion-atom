clear all;

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%%%%%%%%%% Energy ranges in Ba+-Li%%%%%%%%%%%
alpha=164;
alpha_He=1.38;
redmass=1822.8884/(1/138+1/6);
redmassHe=1822.8884/(1/138+1/4);

omega50=1.519e-10*0.05*2*pi;
omega100=1.519e-10*0.1*2*pi;
omega200=1.519e-10*0.2*2*pi;
omega500=1.519e-10*0.5*2*pi;
omega1000=1.519e-10*2*pi;

TrappingPotential=@(m,omega,alpha,l,r) -0.5*alpha*(sqrt(m*omega).*(-2.*exp(-m*omega.*r.^2)./r/sqrt(pi)+erf(sqrt(m*omega)*r)./sqrt(m*omega)./r.^2)).^2+l*(l+1)./2/redmass./r.^2;
Potentialatomion=@(redmass,alpha,l,r) -0.5*alpha./r.^4+l*(l+1)./2/redmass./r.^2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Li_pc=readtable('Li/Li_pt_analytical_and_dif.csv');
omega50SC_Li=readtable('Li/Li_50k_spherical_landau.csv'); 
omega100SC_Li=readtable('Li/Li_100k_spherical_landau.csv');
omega200SC_Li=readtable('Li/Li_200k_spherical_landau.csv');
omega500SC_Li=readtable('Li/Li_500k_spherical_landau.csv');

omega50_Li=readtable('Li/Li_50k_spherical_analytical_and_dif.csv');
omega100_Li=readtable('Li/Li_100k_spherical_analytical_and_dif.csv');
omega200_Li=readtable('Li/Li_200k_spherical_analytical_and_dif.csv');
omega500_Li=readtable('Li/Li_500k_spherical_analytical_and_dif.csv');


He_pc=readtable('He/He_pt_analytical_and_dif.csv');
omega50SC_He=readtable('He/He_50k_spherical_landau.csv');
omega100SC_He=readtable('He/He_100k_spherical_landau.csv');
omega200SC_He=readtable('He/He_200k_spherical_landau.csv');
omega500SC_He=readtable('He/He_500k_spherical_landau.csv');

omega50_He=readtable('He/He_50k_spherical_analytical_and_dif.csv');
omega100_He=readtable('He/He_100k_spherical_analytical_and_dif.csv');
omega200_He=readtable('He/He_200k_spherical_analytical_and_dif.csv');
omega500_He=readtable('He/He_500k_spherical_analytical_and_dif.csv');


Ecol=logspace(-6,2,25)*3.166811e-06;

semiclass_elastic=-4*pi/2/redmass*gamma(-2/3)/3/2^(1/3)*(pi*redmass^2*alpha/4)^(2/3)./Ecol.^(1/3);

semiclass_elastic_He=-4*pi/2/redmassHe*gamma(-2/3)/3/2^(1/3)*(pi*redmassHe^2*alpha_He/4)^(2/3)./Ecol.^(1/3);


figure(10);clf;% panel (b) of Fig 2 of the paper
plot(Li_pc.Var1,Li_pc.Var2,'k','LineWidth',2)
hold on
 plot(omega50_Li.Var1,omega50_Li.Var2,'-r','LineWidth',2)
 hold on
 plot(omega100_Li.Var1,omega100_Li.Var2,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_Li.Var1,omega200_Li.Var2,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_Li.Var1,omega500_Li.Var2,'-b','LineWidth',2)
xlim([1e-6 100])

 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')

 set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel(' $\sigma_e$(a$_0^2$) ','Fontsize',24,'Interpreter','latex')


figure(20);clf;% Fig 5 of the paper

subplot(2,2,1)
plot(Li_pc.Var1,Li_pc.Var3,'k','LineWidth',2)
hold on
 plot(omega50_Li.Var1,omega50_Li.Var3,'-r','LineWidth',2)
 hold on
 plot(omega100_Li.Var1,omega100_Li.Var3,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_Li.Var1,omega200_Li.Var3,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_Li.Var1,omega500_Li.Var3,'-b','LineWidth',2)
xlim([1e-6 100])
 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')

 set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel('$\sigma_d$(a$_0^2$)','Fontsize',24,'Interpreter','latex')

subplot(2,2,3)
plot(Li_pc.Var1,Li_pc.Var3,'k','LineWidth',2)
hold on
 plot(omega50_Li.Var1,omega50_Li.Var4,'-r','LineWidth',2)
 hold on
 plot(omega100_Li.Var1,omega100_Li.Var4,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_Li.Var1,omega200_Li.Var4,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_Li.Var1,omega500_Li.Var4,'-b','LineWidth',2)
 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')
xlim([1e-6 100])
set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel('$\sigma_\eta$(a$_0^2$)','Fontsize',24,'Interpreter','latex')

subplot(2,2,2)
plot(He_pc.Var1,He_pc.Var3,'k','LineWidth',2)
hold on
 plot(omega50_He.Var1,omega50_He.Var3,'-r','LineWidth',2)
 hold on
 plot(omega100_He.Var1,omega100_He.Var3,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_He.Var1,omega200_He.Var3,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_He.Var1,omega500_He.Var3,'-b','LineWidth',2)
xlim([1e-6 100])
 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')

 set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel('$\sigma_d$(a$_0^2$)','Fontsize',24,'Interpreter','latex')

subplot(2,2,4)
plot(He_pc.Var1,He_pc.Var3,'k','LineWidth',2)
hold on
 plot(omega50_He.Var1,omega50_He.Var4,'-r','LineWidth',2)
 hold on
 plot(omega100_He.Var1,omega100_He.Var4,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_He.Var1,omega200_He.Var4,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_He.Var1,omega500_He.Var4,'-b','LineWidth',2)
 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')

 xlim([1e-6 100])
set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel('$\sigma_\eta$(a$_0^2$)','Fontsize',24,'Interpreter','latex')


%%%%%%%%%%%% Average momentum squared transfer %%%%%%%%%%%%%%%

Ecol=omega100_Li.Var1;
redmass=1822.8884/(1/138+1/6);
p=sqrt(2*redmass*Ecol);

omega50=1.519e-10*0.05;
omega100=1.519e-10*0.1;
omega200=1.519e-10*0.2;
omega500=1.519e-10*0.5;
omega1000=1.519e-10*2;


figure(1003);clf; % panel (a) of Fig.6 of the paper
plot(Ecol,(2*Li_pc.Var3)./Li_pc.Var2.*Ecol,'-k','LineWidth',2)
hold on
plot(Ecol,(2*omega50_Li.Var3)./omega50_Li.Var2.*Ecol,'-r','LineWidth',2)
hold on
plot(Ecol,(2*omega100_Li.Var3)./omega100_Li.Var2.*Ecol,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
hold on
plot(Ecol,(2*omega200_Li.Var3)./omega200_Li.Var2.*Ecol,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
hold on
plot(Ecol,(2*omega500_Li.Var3)./omega500_Li.Var2.*Ecol,'-b','LineWidth',2)
     hold on
plot(Ecol,(Ecol).^(5/6)/4,':','LineWidth',2)

 legend('Point Charged','50kHz','100kHz','200kHz','500kHz')
set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
ylabel(' Energy Transfer (K)','Fontsize',24,'Interpreter','latex')
xlabel('$E$(K)','Fontsize',24,'Interpreter','latex')
xlim([1e-6 100])

figure(2001);clf;%Panel (b) of Fig.4 of the paper
plot(He_pc.Var1,He_pc.Var2,'k','LineWidth',2)
hold on
 plot(omega50_He.Var1,omega50_He.Var2,'-r','LineWidth',2)
 hold on
 plot(omega100_He.Var1,omega100_He.Var2,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_He.Var1,omega200_He.Var2,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_He.Var1,omega500_He.Var2,'-b','LineWidth',2)

 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')

 set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel(' $\sigma_e$(a$_0^2$) ','Fontsize',24,'Interpreter','latex')

xlim([1e-6 100])
%%%%%%%%%%%% Average momentum squared transfer %%%%%%%%%%%%%%%

Ecol=omega100_He.Var1;
redmass=1822.8884/(1/138+1/6);
p=sqrt(2*redmass*Ecol);

omega50=1.519e-10*0.05;
omega100=1.519e-10*0.1;
omega200=1.519e-10*0.2;
omega500=1.519e-10*0.5;
omega1000=1.519e-10*2;


figure(2003);clf;%panel (b) of Fig.6 of the paper
plot(Ecol,(2*He_pc.Var3)./He_pc.Var2.*Ecol,'-k','LineWidth',2)
hold on
plot(Ecol,(2*omega50_He.Var3)./omega50_He.Var2.*Ecol,'-r','LineWidth',2)
hold on
plot(Ecol,(2*omega100_He.Var3)./omega100_He.Var2.*Ecol,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
hold on
plot(Ecol,(2*omega200_He.Var3)./omega200_He.Var2.*Ecol,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
hold on
plot(Ecol,(2*omega500_He.Var3)./omega500_He.Var2.*Ecol,'-b','LineWidth',2)
     hold on
plot(Ecol,(Ecol).^(5/6)/2,':','LineWidth',2)
xlim([1e-6 100])
 legend('Point Charged','50kHz','100kHz','200kHz','500kHz')
set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
ylabel(' Energy Transfer Efficiency','Fontsize',24,'Interpreter','latex')
xlabel('$E$(K)','Fontsize',24,'Interpreter','latex')