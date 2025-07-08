fs=20;
ms=8;
dia = 2.9e-3;
aoa = [1e-11, 2e-11, 3e-11, 4e-11, 5e-11, 6e-11, 7e-11, 8e-11];
%%
l0 = 1e-3;
L0 = 0.105;
l = 7;
N = 2^l - 1;
n0 = 1.000275;
Ld = 635e-9;
k0 = 1/Ld * 2 * pi() * n0;
Lx = 0.1;
dx = Lx/(N+1);
delta = dx;
x = ((-(N+1)/2 : (N+1)/2-1)) * delta;
k = linspace(0, 1/dx, 5.*N+1);
turblen = 100;
k=k(2:end);
%Cn2 = linspace(0, 1e-10, 100);
Km = 3.3/l0; % inner scale frequency [1/m]
K0 = 2*pi/L0;
%PSD_phi = 0.033.* exp(-(k./Km).^2)./ (k.^2 + K0^2).^(11/6) .* (1+1.802.*(k./Km) - 0.254 .*(k./Km).^(7/6));
PSD_phi = 0.033.* exp(-(k./Km).^2)./ (k.^2).^(11/6) .* (1+1.802.*(k./Km) - 0.254 .*(k./Km).^(7/6)).*(1-exp(-k.^2/K0.^2));
Km1 = 5.92/l0; % inner scale frequency [1/m]
K01 = 1/L0*2*pi;
PSD_phi1 = 0.033.* exp(-(k./Km1).^2)./ (k.^2 + K01^2).^(11/6);
Km2=5.92/l0;
PSD_phi2 = 0.033*k.^(-11/3).*exp(-(k./Km2).^2);
J = besselj(0,dia.*k);
D_mvk = aoa.*(k0.*dia).^2;
Cn2 = D_mvk ./( 8.*pi^2.*k0.^2*turblen.*sum(k.*PSD_phi.*(1-J),'all').*(k(2)-k(1)) );
Cn2_mvk = D_mvk ./( 8.*pi^2.*k0.^2*turblen.*sum(k.*PSD_phi1.*(1-J),'all').*(k(2)-k(1)) );
Cn2_Tar = D_mvk ./( 8.*pi^2.*k0.^2*turblen.*sum(k.*PSD_phi2.*(1-J),'all').*(k(2)-k(1)) );
Cn2_kol = aoa.*dia.^(1/3)./2.91./turblen;
a= Cn2(end)/aoa(end);
b= Cn2_mvk(end)/aoa(end);
c= Cn2_kol(end)/aoa(end);
d= Cn2_Tar(end)/aoa(end);
shape=["--s","--d","--o"];
shape_modified=["-s","-d","-o"];
shape_exp=["s","d","o"];
figure(1)
load('simulation_result_all_1012_modified_exp_2pi_6mm_l0prop_5mm.mat');
cn2=a/c*2E-13*(5:5:50).^2.0022;
cn2=[1.7970E-14 cn2];
for i=1:3
    plot(cn2,Ecc(i,:),shape(i),'MarkerSize',ms,'LineWidth',1.5);
    hold on;
end 
load("Experimental.mat");
cn2=[1.7970E-14 2E-13*a/c*[15 30 40 50].^2.0022];
for i=1:3
    errorbar(cn2,Ecc_exp(i,:),Ecc_exp_std(i,:),shape_exp(i),'MarkerSize',ms,'MarkerFaceColor','auto','LineWidth',1.5);
    hold on;
end 
hold off;
legend('m=0 modified','m=1 modified','m=2 modified','m=0 exp','m=1 exp','m=2 exp','NumColumns',2,'Location','North');
xlabel('C_{n}^{2} (m^{-2/3})');
ylabel('Eccentricity');
axis([0 0.74E-9 0 1]);
set(gca,'FontSize',fs);
saveas(gca,'ecc_mod_2pi.pdf','pdf');
figure(2)
load('simulation_result_all_1012_modified_exp_2pi_6mm_l0prop_5mm.mat');
cn2=a/c*2E-13*(5:5:50).^2.0022;
cn2=[1.7970E-14 cn2];
for i=1:3
    plot(cn2,Ell(i,:),shape(i),'MarkerSize',ms,'LineWidth',1.5);
    hold on;
end 
load("Experimental.mat");
cn2=[1.7970E-14 2E-13*a/c*[15 30 40 50].^2.0022];
for i=1:3
    errorbar(cn2,Ell_exp(i,:),Ell_exp_std(i,:),shape_exp(i),'MarkerSize',ms,'MarkerFaceColor','auto','LineWidth',1.5);
    hold on;
end 
hold off;
legend('m=0 modified','m=1 modified','m=2 modified','m=0 exp','m=1 exp','m=2 exp','NumColumns',2,'Location','South');
xlabel('C_{n}^{2} (m^{-2/3})');
ylabel('Ellipticity');
axis([0 0.74E-9 0.6 1]);
set(gca,'FontSize',fs);
saveas(gca,'ell_mod_2pi.pdf','pdf');