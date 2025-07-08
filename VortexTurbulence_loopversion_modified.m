%% Initialize Matlab
clc
clear all
% clear variables
close all
dia = 2.9e-3;
aoa = [1e-11, 2e-11, 3e-11, 4e-11, 5e-11, 6e-11, 7e-11, 8e-11];
l0 = 1e-3;
L0 = 0.105;
l0_prop=5E-3;
L0_prop=1;
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
turblen = L0/0.7;
k=k(2:end);
%Cn2 = linspace(0, 1e-10, 100);
Km = 3.3/l0; % inner scale frequency [1/m]
K0 = 8*pi/L0;
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
Km = 3.3/l0_prop; % inner scale frequency [1/m]
K0 = 8*pi/L0_prop;
%PSD_phi = 0.033.* exp(-(k./Km).^2)./ (k.^2 + K0^2).^(11/6) .* (1+1.802.*(k./Km) - 0.254 .*(k./Km).^(7/6));
PSD_phi = 0.033.* exp(-(k./Km).^2)./ (k.^2).^(11/6) .* (1+1.802.*(k./Km) - 0.254 .*(k./Km).^(7/6)).*(1-exp(-k.^2/K0.^2));
Km1 = 5.92/l0_prop; % inner scale frequency [1/m]
K01 = 1/L0_prop*2*pi;
PSD_phi1 = 0.033.* exp(-(k./Km1).^2)./ (k.^2 + K01^2).^(11/6);
Km2=5.92/l0_prop;
PSD_phi2 = 0.033*k.^(-11/3).*exp(-(k./Km2).^2);
J = besselj(0,dia.*k);
D_mvk = aoa.*(k0.*dia).^2;
Cn2 = D_mvk ./( 8.*pi^2.*k0.^2*turblen.*sum(k.*PSD_phi.*(1-J),'all').*(k(2)-k(1)) );
Cn2_mvk = D_mvk ./( 8.*pi^2.*k0.^2*turblen.*sum(k.*PSD_phi1.*(1-J),'all').*(k(2)-k(1)) );
Cn2_Tar = D_mvk ./( 8.*pi^2.*k0.^2*turblen.*sum(k.*PSD_phi2.*(1-J),'all').*(k(2)-k(1)) );
Cn2_kol = aoa.*dia.^(1/3)./2.91./turblen;
a1= Cn2(end)/aoa(end);
b1= Cn2_mvk(end)/aoa(end);
c1= Cn2_kol(end)/aoa(end);
d1= Cn2_Tar(end)/aoa(end);
%dt=[15 30 40 50];
dt=5:5:50;
cn2=2E-13*dt.^2.0022;
cn2=cn2*a/c;
prop=1E-13*a1/c1;
cn2=[prop cn2];
lambda=532E-9;
k=2*pi/lambda;
Dz_turb=0.15;
Dz=100;
turb_nscr=2;
prop_nscr=5;
R0=zeros(turb_nscr,length(cn2));
for i=1:length(cn2)
    R0(:,i)=Cn2r0(Dz_turb,k,cn2(i),turb_nscr);
end
R0_prop=Cn2r0(Dz,k,cn2(1),prop_nscr);
Ecc=zeros(3,length(R0));Eccstd=zeros(3,length(R0));
Ell=zeros(3,length(R0));Ellstd=zeros(3,length(R0));
SI=zeros(3,length(R0));
incrementNumber = 2;
incrementNumber_prop = 5;
realizationNumber = 1; % how many frames of realizations we could have
frameNumber=100;
Eccentricity = zeros(1,frameNumber);
Ellipticity = zeros(1,frameNumber);
Imax=zeros(1,frameNumber);      
    %% Parameters of field griding and propagation
    
    k0 = 1i * 2 * pi / lambda;    % complex wave vector
    scale_factor = 1.0;
    w0 = 0.3E-2; % beam waist radius
    delta = 200E-6 * scale_factor;
    N=1499;
    wmax=delta*N;
    grid1D_point = round(wmax/delta+1);% the amount of points of linear dimension                            
    warning('wmax and delta') % Anywhere Door                        
      
    [X,Y] = meshgrid(-wmax/2:delta:wmax/2,-wmax/2:delta:wmax/2);
    
    [Theta,Rho] = cart2pol(X,Y);% transfer from Cartesian coordinate to polar coordinate
                                    
    x = X(1,:);                 % take 1D vector from the 2D grid, row
    y = Y(:,1);                 % take 1D vector from the 2D grid, column
    %% Parameter of k space
    
    z = 100; 
    
    %% Parameter of Frequency filter (analytical expression), for propagation purpose
    
    WMAX = grid1D_point/wmax;	% length edge of the light field in spatial 
                                % frequency domain
                                % with the same number of points
    
    DELTA = 1/wmax;             % step size in spatial frequency domain
    
    [fX,fY] = meshgrid...
        ((-WMAX/2):DELTA:(WMAX/2-DELTA),(-WMAX/2):DELTA:(WMAX/2-DELTA));
                                % regrid the spatial frequency domain
for i=0:1:2
    for j=1:length(cn2)
        disp("looping number: "+num2str(i*length(R0)+j))
        for frameindex=1:frameNumber   
            disp(frameindex)
        chargeNumber = i;
        phaseVortex = (Theta) * chargeNumber;
        U_Gauss_x = exp(-((0.916*X).^2/(w0^2)+Y.^2/(w0^2)))...
            .*exp(1j * phaseVortex);
        N = grid1D_point; % number of grid points per side
        L0 = 0.105; % outer scale [m]
        l0 = 0.001;% inner scale [m]
        delta = D/N; % grid spacing [m]

        %% Turbulence in the beginning
        % zi = 20* meters/incrementNumber;
        zi = Dz_turb/incrementNumber;
        zi_prop=z/incrementNumber_prop;
        ZVAL_RS_zi = exp(k0*zi*sqrt(1-(lambda*fX).^2-(lambda*fY).^2)); 
        ZVAL_RS_prop = exp(k0*zi_prop*sqrt(1-(lambda*fX).^2-(lambda*fY).^2));

        for i_realization_index = 1 : realizationNumber
            for i_index = 1 : incrementNumber

                r0=R0(i_index,j);
                [phz_lo, phz_hi] = ft_sh_phase_screen_modified_exp(r0, N, delta, L0, l0);
                phz = phz_lo + phz_hi;
                phase_screen = exp(1j*mod(phz,2*pi));
                U_Gauss_zi_x= ifft2...
                    ((fft2(U_Gauss_zi_x .* phase_screen)/grid1D_point).*...
                    fftshift(ZVAL_RS_zi)).*grid1D_point; 
            end
        
            for i_index = 1 : incrementNumber_prop
                    r0_prop=R0_prop(i_index);
                    [phz_lo, phz_hi] = ft_sh_phase_screen_modified_exp(r0_prop, N, delta, L0_prop, l0_prop);
                    phz = phz_lo + phz_hi;
                phase_screen = exp(1j*mod(phz,2*pi));
                U_Gauss_zi_x= ifft2...
                    ((fft2(U_Gauss_zi_x .* phase_screen)/grid1D_point).*...
                    fftshift(ZVAL_RS_prop)).*grid1D_point;
            end
            
            I_Gauss_zi_realization = abs(U_Gauss_zi_x).^2;
        end
        bw = zeros(grid1D_point,grid1D_point);
        for i_index = 1:realizationNumber
            bw=I_Gauss_zi_realization;
            totalx=0;totaly=0;totalIxx=0;totalIyy=0;totalIxy=0;
            for p=1:N
                for q=1:N
                    totalx=totalx+p*bw(p,q);
                    totaly=totaly+q*bw(p,q);
                end
            end
            totalI=sum(sum(bw));
            cx=round(totalx/totalI);cy=round(totaly/totalI);
            for p=1:N
                for q=1:N
                    totalIxx=totalIxx+(p-cx)^2*bw(p,q);
                    totalIyy=totalIyy+(q-cy)^2*bw(p,q);
                    totalIxy=totalIxy+(p-cx)*(q-cy)*bw(p,q);
                end
            end
            muxx=totalIxx/totalI;muyy=totalIyy/totalI;muxy=totalIxy/totalI;
            orientation=1/2*atan((2*muxy)/(muxx-muyy)); %orientation in rad
            dd=sqrt(4*muxy^2+(muxx-muyy)^2);
            lambda1=((muxx+muyy)+dd)/2;lambda2=((muxx+muyy)-dd)/2;
            sigmaM=4*sqrt(abs(lambda1));
            sigmam=4*sqrt(abs(lambda2));
            Eccentricity(frameindex)=sqrt(sigmaM^2-sigmam^2)/sigmaM;
            Ellipticity(frameindex)=sigmam/sigmaM;
            r=200;
            IAzi=zeros(1,r+1);
            deltaa=0.2;stepp=0.2;
            for ii=0:r
                for s=ii-deltaa:stepp:ii+deltaa
                    for theta=1:360
                        xx=round(cy+s*cos(theta*pi/180));
                        yy=round(cx+s*sin(theta*pi/180));
                        IAzi(ii+1)=IAzi(ii+1)+bw(xx,yy);
                    end
                end
                IAzi(ii+1)=IAzi(ii+1)/(2*deltaa/stepp+1);
            end
            Imax(frameindex)=max(IAzi);
        end
        end
        ecc=mean(Eccentricity);
        eccstd=std(Eccentricity);
        ell=mean(Ellipticity);
        ellstd=std(Ellipticity);
        Ecc(i+1,j)=ecc;Eccstd(i+1,j)=eccstd;
        Ell(i+1,j)=ell;Ellstd(i+1,j)=ellstd;
        SI(i+1,j)=(mean(Imax.^2)-mean(Imax)^2)/(mean(Imax)^2);
    end
end

%% Plot
shape=["s","d","o"];
figure(1)
for i=1:3
    errorbar(cn2,Ecc(i,:),Eccstd(i,:),shape(i),'MarkerSize',8,'LineWidth',1.5);
    hold on;
end 
hold off;
legend('m=0','m=1','m=2','Location','Best');
xlabel('C_{n}^{2} (m^{-2/3})');
ylabel('Eccentricity');
axis([0 1.1E-9 0 0.8]);
set(gca,'FontSize',16);
set(gca,'FontWeight','bold');
figure(2)
for i=1:3
    errorbar(cn2,Ell(i,:),Ellstd(i,:),shape(i),'MarkerSize',8,'LineWidth',1.5);
    hold on;
end 
hold off;
legend('m=0','m=1','m=2','Location','Best');
xlabel('C_{n}^{2} (m^{-2/3})');
ylabel('Ellipticity');
axis([0 1.1E-9 0.6 1]);
set(gca,'FontSize',16);
set(gca,'FontWeight','bold');
figure(3)
color=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];
for i=1:3
    plot(cn2,SI(i,:),'o','MarkerSize',10,'MarkerFaceColor',color(i,:),'MarkerEdgeColor',color(i,:));hold on;
end
hold off;
legend('m=0','m=1','m=2','Location','Best');
xlabel('C_n^2 (m^{-2/3})');
ylabel('Scintillation Index');
set(gca,'FontSize',16);
set(gca,'FontWeight','bold');
axis([0 1.1E-9 0 0.64]);
save('simulation_result_all_0920_modified_exp_8pi_6mm_l0prop_5mm.mat','Ecc','Eccstd','Ell','Ellstd','SI','cn2');