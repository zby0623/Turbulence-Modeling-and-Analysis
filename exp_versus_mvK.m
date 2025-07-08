fs=20;
ms=8;
shape=["--s","--d","--o"];
shape_modified=["-s","-d","-o"];
shape_exp=["s","d","o"];
figure(1)
load("simulation_result_all_1009_mvK.mat");
cn2=0.6752/0.3267*2E-13*(5:5:50).^2.0022;
cn2=[1.7970E-14 cn2];
for i=1:3
    plot(cn2,Ecc(i,:),shape(i),'MarkerSize',ms,'LineWidth',1.5);
    hold on;
end 
load("Experimental.mat");
cn2=[1.7970E-14 0.6752/0.3267*2E-13*[15 30 40 50].^2.0022];
for i=1:3
    errorbar(cn2,Ecc_exp(i,:),Ecc_exp_std(i,:),shape_exp(i),'MarkerSize',ms,'MarkerFaceColor','auto','LineWidth',1.5);
    hold on;
end 
hold off;
legend('m=0 mvK','m=1 mvK','m=2 mvK','m=0 exp','m=1 exp','m=2 exp','NumColumns',2,'Location','North');
xlabel('C_{n}^{2} (m^{-2/3})');
ylabel('Eccentricity');
axis([0 1.05E-9 0 1]);
set(gca,'FontSize',fs);
saveas(gca,'ecc_mvK.pdf','pdf');
figure(2)
load("simulation_result_all_1009_mvK.mat");
cn2=0.6752/0.3267*2E-13*(5:5:50).^2.0022;
cn2=[1.7970E-14 cn2];
for i=1:3
    plot(cn2,Ell(i,:),shape(i),'MarkerSize',ms,'LineWidth',1.5);
    hold on;
end 
load("Experimental.mat");
cn2=[1.7970E-14 0.6752/0.3267*2E-13*[15 30 40 50].^2.0022];
for i=1:3
    errorbar(cn2,Ell_exp(i,:),Ell_exp_std(i,:),shape_exp(i),'MarkerSize',ms,'MarkerFaceColor','auto','LineWidth',1.5);
    hold on;
end 
hold off;
legend('m=0 mvK','m=1 mvK','m=2 mvK','m=0 exp','m=1 exp','m=2 exp','NumColumns',2,'Location','South');
xlabel('C_{n}^{2} (m^{-2/3})');
ylabel('Ellipticity');
axis([0 1.05E-9 0.6 1]);
set(gca,'FontSize',fs);
saveas(gca,'ell_mvK.pdf','pdf');