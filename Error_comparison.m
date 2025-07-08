fs=20;
sz=9;
figure(1)
load('Experimental.mat');
load('simulation_result_all_1009_modified_exp_8pi.mat');
Err_ell_mod=abs(Ell(:,[4 7 9 11])-Ell_exp(:,2:end));
hold on;
load('simulation_result_all_1009_mvK.mat');
Err_ell_mvK=abs(Ell(:,[4 7 9 11])-Ell_exp(:,2:end));
hold on;
load('simulation_result_all_1009_Tatarskii.mat');
Err_ell_tar=abs(Ell(:,[4 7 9 11])-Ell_exp(:,2:end));
hold on;
load('simulation_result_all_1009_Kolmogorov.mat');
Err_ell_kol=abs(Ell(:,[4 7 9 11])-Ell_exp(:,2:end));
hold off;
shape=['s','d','o'];
for i=1:3
    plot(Err_ell_mod(i,:),'LineWidth',1.5,'Marker',shape(i),'Color',"#FF2908",'MarkerSize',sz);hold on;
    plot(Err_ell_mvK(i,:),'LineWidth',1.5,'Marker',shape(i),'Color',"#FCA800",'MarkerSize',sz);hold on;
    plot(Err_ell_tar(i,:),'LineWidth',1.5,'Marker',shape(i),'Color',"#1D00DB",'MarkerSize',sz);hold on;
    plot(Err_ell_kol(i,:),'LineWidth',1.5,'Marker',shape(i),'Color',"#00A204",'MarkerSize',sz);hold on;
end
hold off;
set(gca,'FontSize',fs);
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',[15 30 40 50]);
xlabel('Temperature Difference (^{o}C)');
ylabel('Ellipticity Error');
saveas(gca,'err_ell.pdf','pdf');
figure(2)
load('Experimental.mat');
load('simulation_result_all_1009_modified_exp_8pi.mat');
Err_ecc_mod=abs(Ecc(:,[4 7 9 11])-Ecc_exp(:,2:end));
hold on;
load('simulation_result_all_1009_mvK.mat');
Err_ecc_mvK=abs(Ecc(:,[4 7 9 11])-Ecc_exp(:,2:end));
hold on;
load('simulation_result_all_1009_Tatarskii.mat');
Err_ecc_tar=abs(Ecc(:,[4 7 9 11])-Ecc_exp(:,2:end));
hold on;
load('simulation_result_all_1009_Kolmogorov.mat');
Err_ecc_kol=abs(Ecc(:,[4 7 9 11])-Ecc_exp(:,2:end));
hold off;
shape=['s','d','o'];
for i=1:3
    plot(Err_ecc_mod(i,:),'LineWidth',1.5,'Marker',shape(i),'Color',"#FF2908",'MarkerSize',sz);hold on;
    plot(Err_ecc_mvK(i,:),'LineWidth',1.5,'Marker',shape(i),'Color',"#FCA800",'MarkerSize',sz);hold on;
    plot(Err_ecc_tar(i,:),'LineWidth',1.5,'Marker',shape(i),'Color',"#1D00DB",'MarkerSize',sz);hold on;
    plot(Err_ecc_kol(i,:),'LineWidth',1.5,'Marker',shape(i),'Color',"#00A204",'MarkerSize',sz);hold on;
end
hold off;
set(gca,'FontSize',fs);
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',[15 30 40 50]);
xlabel('Temperature Difference (^{o}C)');
ylabel('Eccentricity Error');
saveas(gca,'err_ecc.pdf','pdf');
mean_ell=zeros(1,4);
mean_ell(1)=mean(mean(Err_ell_mod));
mean_ell(2)=mean(mean(Err_ell_mvK));
mean_ell(3)=mean(mean(Err_ell_tar));
mean_ell(4)=mean(mean(Err_ell_kol));
mean_ecc=zeros(1,4);
mean_ecc(1)=mean(mean(Err_ecc_mod));
mean_ecc(2)=mean(mean(Err_ecc_mvK));
mean_ecc(3)=mean(mean(Err_ecc_tar));
mean_ecc(4)=mean(mean(Err_ecc_kol));



