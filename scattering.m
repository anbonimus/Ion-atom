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

Li_pc=readtable('C:\Users\12107\OneDrive\Desktop\Li\Li_point_charge_cross_section.csv');
omega50_Li=readtable('C:\Users\12107\OneDrive\Desktop\Li\Li_50kHz_cross_section.csv');
omega100_Li=readtable('C:\Users\12107\OneDrive\Desktop\Li\Li_100kHz_cross_section.csv');
omega200_Li=readtable('C:\Users\12107\OneDrive\Desktop\Li\Li_200kHz_cross_section.csv');
omega500_Li=readtable('C:\Users\12107\OneDrive\Desktop\Li\Li_500kHz_cross_section.csv');
omega1000_Li=readtable('C:\Users\12107\OneDrive\Desktop\Li\Li_1MHz_cross_section.csv');

He_pc=readtable('He/He_point_charge_cross_section.csv');
omega50_He=readtable('He/He_50kHz_cross_section.csv');
omega100_He=readtable('He/He_100kHz_cross_section.csv');
omega200_He=readtable('He/He_200kHz_cross_section.csv');
omega500_He=readtable('He/He_500kHz_cross_section.csv');
omega1000_He=readtable('He/He_1MHz_cross_section.csv');

Ecol=logspace(-6,2,25)*3.166811e-06;

semiclass_elastic=-2*pi/redmass*gamma(-2/3)/3/2^(1/3)*(pi*redmass^2*alpha/4)^(2/3)./Ecol.^(1/3);

semiclass_elastic_He=-2*pi/redmassHe*gamma(-2/3)/3/2^(1/3)*(pi*redmassHe^2*alpha_He/4)^(2/3)./Ecol.^(1/3);


figure(10);clf;% panel (b) of Fig 2 of the paper
plot(Li_pc.r_a0_,Li_pc.elastic,'k','LineWidth',2)
hold on
 plot(omega50_Li.r_a0_,omega50_Li.elastic,'-r','LineWidth',2)
 hold on
 plot(omega100_Li.r_a0_,omega100_Li.elastic,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_Li.r_a0_,omega200_Li.elastic,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_Li.r_a0_,omega500_Li.elastic,'-b','LineWidth',2)
xlim([1e-6 100])

 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')

 set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel(' $\sigma_e$(a$_0^2$) ','Fontsize',24,'Interpreter','latex')

figure(10);clf;% panel (b) of Fig 2 of the paper
plot(He_pc.r_a0_,He_pc.elastic,'k','LineWidth',2)
hold on
 plot(omega50_He.r_a0_,omega50_He.elastic,'-r','LineWidth',2)
 hold on
 plot(omega100_He.r_a0_,omega100_He.elastic,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_He.r_a0_,omega200_He.elastic,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_He.r_a0_,omega500_He.elastic,'-b','LineWidth',2)
xlim([1e-6 100])

 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')

 set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel(' $\sigma_e$(a$_0^2$) ','Fontsize',24,'Interpreter','latex')



figure(20);clf;% Fig 5 of the paper

subplot(2,2,1)
plot(Li_pc.r_a0_,Li_pc.diffusion,'k','LineWidth',2)
hold on
 plot(omega50_Li.r_a0_,omega50_Li.diffusion,'-r','LineWidth',2)
 hold on
 plot(omega100_Li.r_a0_,omega100_Li.diffusion,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_Li.r_a0_,omega200_Li.diffusion,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_Li.r_a0_,omega500_Li.diffusion,'-b','LineWidth',2)
xlim([1e-6 100])
 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')

 set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel('$\sigma_d$(a$_0^2$)','Fontsize',24,'Interpreter','latex')

subplot(2,2,3)
plot(Li_pc.r_a0_,Li_pc.viscosity,'k','LineWidth',2)
hold on
 plot(omega50_Li.r_a0_,omega50_Li.viscosity,'-r','LineWidth',2)
 hold on
 plot(omega100_Li.r_a0_,omega100_Li.viscosity,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_Li.r_a0_,omega200_Li.viscosity,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_Li.r_a0_,omega500_Li.viscosity,'-b','LineWidth',2)
 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')
xlim([1e-6 100])
set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel('$\sigma_\eta$(a$_0^2$)','Fontsize',24,'Interpreter','latex')

subplot(2,2,2)
plot(He_pc.r_a0_,He_pc.diffusion,'k','LineWidth',2)
hold on
 plot(omega50_He.r_a0_,omega50_He.diffusion,'-r','LineWidth',2)
 hold on
 plot(omega100_He.r_a0_,omega100_He.diffusion,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_He.r_a0_,omega200_He.diffusion,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_He.r_a0_,omega500_He.diffusion,'-b','LineWidth',2)
xlim([1e-6 100])
 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')

 set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel('$\sigma_d$(a$_0^2$)','Fontsize',24,'Interpreter','latex')

subplot(2,2,4)
plot(He_pc.r_a0_,He_pc.viscosity,'k','LineWidth',2)
hold on
 plot(omega50_He.r_a0_,omega50_He.viscosity,'-r','LineWidth',2)
 hold on
 plot(omega100_He.r_a0_,omega100_He.viscosity,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
 hold on
 plot(omega200_He.r_a0_,omega200_He.viscosity,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
 hold on
 plot(omega500_He.r_a0_,omega500_He.viscosity,'-b','LineWidth',2)
 legend('Point Charge','50kHz','100kHz','200kHz','500kHz')

 xlim([1e-6 100])
set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
xlabel('$E$ (K)','Fontsize',24,'Interpreter','latex')
ylabel('$\sigma_\eta$(a$_0^2$)','Fontsize',24,'Interpreter','latex')



%%%%%%%%%%%% Average momentum squared transfer %%%%%%%%%%%%%%%

Ecol=omega100_Li.r_a0_;
redmass=1822.8884/(1/138+1/6);
p=sqrt(2*redmass*Ecol);

omega50=1.519e-10*0.05;
omega100=1.519e-10*0.1;
omega200=1.519e-10*0.2;
omega500=1.519e-10*0.5;
omega1000=1.519e-10*2;


figure(1003);clf; % panel (a) of Fig.6 of the paper
plot(Ecol,(2*Li_pc.diffusion)./Li_pc.elastic.*Ecol,'-k','LineWidth',2)
hold on
plot(Ecol,(2*omega50_Li.diffusion)./omega50_Li.elastic.*Ecol,'-r','LineWidth',2)
hold on
plot(Ecol,(2*omega100_Li.diffusion)./omega100_Li.elastic.*Ecol,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
hold on
plot(Ecol,(2*omega200_Li.diffusion)./omega200_Li.elastic.*Ecol,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
hold on
plot(Ecol,(2*omega500_Li.diffusion)./omega500_Li.elastic.*Ecol,'-b','LineWidth',2)
     hold on
plot(Ecol,(Ecol).^(5/6)/4,':','LineWidth',2)

 legend('Point Charged','50kHz','100kHz','200kHz','500kHz')
set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
ylabel(' Energy Transfer (K)','Fontsize',24,'Interpreter','latex')
xlabel('$E$(K)','Fontsize',24,'Interpreter','latex')
xlim([1e-6 100])

%%%%%%%%%%%% Average momentum squared transfer %%%%%%%%%%%%%%%

Ecol=omega100_He.r_a0_;
redmass=1822.8884/(1/138+1/6);
p=sqrt(2*redmass*Ecol);

omega50=1.519e-10*0.05;
omega100=1.519e-10*0.1;
omega200=1.519e-10*0.2;
omega500=1.519e-10*0.5;
omega1000=1.519e-10*2;


figure(2003);clf;%panel (b) of Fig.6 of the paper
plot(Ecol,(2*He_pc.diffusion)./He_pc.elastic.*Ecol,'-k','LineWidth',2)
hold on
plot(Ecol,(2*omega50_He.diffusion)./omega50_He.elastic.*Ecol,'-r','LineWidth',2)
hold on
plot(Ecol,(2*omega100_He.diffusion)./omega100_He.elastic.*Ecol,'-g','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
hold on
plot(Ecol,(2*omega200_He.diffusion)./omega200_He.elastic.*Ecol,'-m','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
hold on
plot(Ecol,(2*omega500_He.diffusion)./omega500_He.elastic.*Ecol,'-b','LineWidth',2)
     hold on
plot(Ecol,(Ecol).^(5/6)/2,':','LineWidth',2)
xlim([1e-6 100])
 legend('Point Charged','50kHz','100kHz','200kHz','500kHz')
set( gca,'FontSize',24,'FontName', 'Helvetica','YScale', 'log','XScale', 'log'  );
ylabel(' Energy Transfer Efficiency','Fontsize',24,'Interpreter','latex')
xlabel('$E$(K)','Fontsize',24,'Interpreter','latex')
