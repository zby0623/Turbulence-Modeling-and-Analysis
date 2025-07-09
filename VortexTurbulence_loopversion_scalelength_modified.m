%% Initialize Matlab
clc
clear 
% clear variables
close all
cn2=2E-10;
lambda=532E-9;
k=2*pi/lambda;
prop=7.19E-15;
Dz_turb=0.15;
Dz=100;
turb_nscr=2;
prop_nscr=5;
L0_prop=1;
l0_prop=0.005;
L0_t=0.1:0.1:1;
l0_t=0.001:0.001:0.01;
R0=zeros(turb_nscr,length(cn2));
for i=1:length(cn2)
    R0(:,i)=Cn2r0(Dz_turb,k,cn2(i),turb_nscr);
end
R0_prop=Cn2r0(Dz,k,prop,prop_nscr);
Ecc_L0=zeros(3,length(L0_t));Eccstd_L0=zeros(3,length(L0_t));
Ell_L0=zeros(3,length(L0_t));Ellstd_L0=zeros(3,length(L0_t));
SI_L0=zeros(3,length(L0_t));
Ecc_l0=zeros(3,length(L0_t));Eccstd_l0=zeros(3,length(L0_t));
Ell_l0=zeros(3,length(L0_t));Ellstd_l0=zeros(3,length(L0_t));
SI_l0=zeros(3,length(L0_t));
incrementNumber = 2;
incrementNumber_prop = 5;
realizationNumber = 1; % how many frames of realizations we could have
frameNumber=200;
Eccentricity_L0 = zeros(1,frameNumber);
Ellipticity_L0 = zeros(1,frameNumber);
Eccentricity_l0 = zeros(1,frameNumber);
Ellipticity_l0 = zeros(1,frameNumber);
Imax_L0=zeros(1,frameNumber);
Imax_l0=zeros(1,frameNumber);
%% Assign unit
microns = 1; % determine scale of units
% microns = 1e-6; % determine scale of units; phase plate generation?

global nanometers nm micrometers um millimetre mm cm meters...
    degrees deg rad inch;
% wUnits(microns)     % use function of zunits to assign units
nanometers	= 1e-3*microns;
nm          = nanometers;
micrometers = microns;
um          = micrometers;
millimetre  = 1e3*microns;
mm          = millimetre;
cm          = 1e4*microns;
meters      = 1e6*microns; % m is the topological charge number
degrees     = atan(microns/microns)/45;
deg         = degrees;
rad         = (atan(microns/microns)/45)*180/pi;
inch        = 25400 * microns;

%% Welcome Interface
fprintf(['Welcom to the program of \n',...
    'VortexTurbulence_20220420.m\n',...
    '',...
    '''*'' means default setting \n']);
    % print welcome words on command window, hit of default prameter
    % setting

lambda = 532 * nm;  % central wavelegth

%% Parameters of field griding and propagation

k0 = 1i * 2 * pi / lambda;    % complex wave vector
scale_factor = 1.0;
w0 = 0.175 *cm; % beam waist radius
delta = 100*um * scale_factor;
N=2047;
wmax=delta*N;
grid1D_point = round(wmax/delta+1);% the amount of points of linear dimension                            
warning('wmax and delta') % Anywhere Door                        
  
[X,Y] = meshgrid(-wmax/2:delta:wmax/2,-wmax/2:delta:wmax/2);
[Theta,Rho] = cart2pol(X,Y);% transfer from Cartesian coordinate to  

%% Parameter of k space

z = 100* meters;


%% Parameter of Frequency filter (analytical expression), for propagation purpose

WMAX = grid1D_point/wmax;	% length edge of the light field in spatial 
                            % frequency domain
                            % with the same number of points

DELTA = 1/wmax;             % step size in spatial frequency domain

[fX,fY] = meshgrid...
    ((-WMAX/2):DELTA:(WMAX/2-DELTA),(-WMAX/2):DELTA:(WMAX/2-DELTA));
                            % regrid the spatial frequency domain
%% Scale bar initialization
scaleBar    = zeros(grid1D_point,grid1D_point);
scaleBar(grid1D_point/2 + round(600*um/delta) : ...
         grid1D_point/2 + round(600*um/delta) + 2,...
         grid1D_point/2 + floor(400*um/delta) : ...
         grid1D_point/2 + ceil (600*um/delta))          = 1;  
for i=0:1:2
    for j=1:length(L0_t)
        disp("looping number: "+num2str(i*length(cn2)+j))
        for frameindex=1:frameNumber
            disp(frameindex)
            chargeNumber = i;
            phaseVortex = (Theta) * chargeNumber;
            pupil=(sqrt(2)*sqrt((X).^2+Y.^2)/w0).^i;
            % phaseVortex = Theta * 2; % charge 2
            U_Gauss_x =exp(-((X).^2/(w0^2)+Y.^2/(w0^2)))...
            .*exp(1j * phaseVortex); 
            %U_Gauss_x=U_Gauss_x.*pupil;
            
            I_Gauss = abs(U_Gauss_x.^2);
            P_Gauss = sum(I_Gauss(:)); % total power of the beam
            I_Gauss_weight = I_Gauss/P_Gauss;
            phase_Gauss = angle(U_Gauss_x);
            
            % figure;imagesc(x/mm,y/mm,I_Gauss);axis image;xlabel('x/mm');ylabel('y/mm');
            %     caxis([0 1.0*max(I_Gauss(:))]); % 1.5
            %     set(gca,'fontSize',15);
            %     Showlimit = 50;xlim([-Showlimit Showlimit]);ylim([-Showlimit Showlimit])% 100
            
            % figure;imagesc(x/mm,y/mm,angle(U_Gauss(:,:,1)));axis image;xlabel('x/mm');ylabel('y/mm')
            %     set(gca,'fontSize',15);
            %     Showlimit = 50;xlim([-Showlimit Showlimit]);ylim([-Showlimit Showlimit])% 100
            
            %% Propagagte Vortex Gaussian
            % z = 20* meters;
            
            ZVAL_RS_z = exp(k0*z*sqrt(1-(lambda*fX).^2-(lambda*fY).^2)); 
                                        % analytical expression of the Fourier
                                        % transform of the Rayleigh-Sommerfeld kernel
            
            % U_Gauss_z = ifftshift(ifft2...
            %     (ifftshift(fftshift(fft2(fftshift((U_Gauss)/grid1D_point).*fftshift(ZVAL_RS_z)))))).*grid1D_point;                             
            U_Gauss_z_x = ifft2...
                ((fft2(U_Gauss_x)/grid1D_point).*fftshift(ZVAL_RS_z)).*grid1D_point; 
            I_Gauss_z = abs(U_Gauss_z_x.^2);
            
            D = wmax/meters; % length of one side of square phase screen [m]
            N = grid1D_point; % number of grid points per side

            %% loop of the inner scale
            L0 = L0_t(1); % outer scale [m]
            l0 = l0_t(j);% inner scale [m]
            delta = D/N; % grid spacing [m]
            %% Turbulence in the beginning
            zi = 0.15*meters/incrementNumber;
            zi_prop=z/incrementNumber_prop;
            ZVAL_RS_zi = exp(k0*zi*sqrt(1-(lambda*fX).^2-(lambda*fY).^2)); 
            ZVAL_RS_prop = exp(k0*zi_prop*sqrt(1-(lambda*fX).^2-(lambda*fY).^2));
                                        % analytical expression of the Fourier
                                        % transform of the Rayleigh-Sommerfeld kernel
            U_Gauss_zi = zeros(grid1D_point,grid1D_point,incrementNumber+incrementNumber_prop+1);  
            U_Gauss_zi_x(:,:,1)= U_Gauss_x; 
            U_Gauss_zi_realization_x = zeros(grid1D_point,grid1D_point,realizationNumber);
            I_Gauss_zi_realization = zeros(grid1D_point,grid1D_point,realizationNumber);
            
            for i_realization_index = 1 : realizationNumber
                for i_index = 1 : incrementNumber
                    r0=R0(i_index);
                    [phz_lo, phz_hi] = ft_sh_phase_screen_modified_exp(r0, N, delta, L0, l0);
                        phz = phz_lo + phz_hi;
                    phase_screen = exp(1j*mod(phz,2*pi));
            
                    U_Gauss_zi_x(:,:,i_index+1)= ifft2...
                        ((fft2(U_Gauss_zi_x(:,:,i_index) .* phase_screen)/grid1D_point).*...
                        fftshift(ZVAL_RS_zi)).*grid1D_point;     
            
                end
            
                for i_index = 1 : incrementNumber_prop
                    r0_prop=R0_prop(i_index);
                        [phz_lo, phz_hi] = ft_sh_phase_screen_modified_exp(r0_prop, N, delta, L0_prop, l0_prop);
                        phz = phz_lo + phz_hi;
                    phase_screen = exp(1j*mod(phz,2*pi));
                
                    U_Gauss_zi_x(:,:,incrementNumber+i_index+1)= ifft2...
                        ((fft2(U_Gauss_zi_x(:,:,incrementNumber+i_index) .* phase_screen)/grid1D_point).*...
                        fftshift(ZVAL_RS_prop)).*grid1D_point;  
                end
            
                U_Gauss_zi_realization_x(:,:,i_realization_index) = U_Gauss_zi_x(:,:,end);
                I_Gauss_zi_realization(:,:,i_realization_index) = abs(U_Gauss_zi_realization_x(:,:,i_realization_index).^2);
            end
            
            I_total = zeros(grid1D_point,grid1D_point);
            for i_index = 1:1:realizationNumber
                I_total = I_total + I_Gauss_zi_realization(:,:,i_index);
            end
            bw = zeros(grid1D_point,grid1D_point);
            for i_index = 1:realizationNumber
                bw=I_Gauss_zi_realization(:,:,i_index);
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
                Eccentricity_l0(frameindex)=sqrt(sigmaM^2-sigmam^2)/sigmaM;
                Ellipticity_l0(frameindex)=sigmam/sigmaM;
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
                Imax_l0(frameindex)=max(IAzi);
            end
        
        end
        ecc_L0=mean(Eccentricity_L0);
        eccstd_L0=std(Eccentricity_L0);
        ell_L0=mean(Ellipticity_L0);
        ellstd_L0=std(Ellipticity_L0);
        Ecc_L0(i+1,j)=ecc_L0;Eccstd_L0(i+1,j)=eccstd_L0;
        Ell_L0(i+1,j)=ell_L0;Ellstd_L0(i+1,j)=ellstd_L0;
        SI_L0(i+1,j)=(mean(Imax_L0.^2)-mean(Imax_L0)^2)/(mean(Imax_L0)^2);
        ecc_l0=mean(Eccentricity_l0);
        eccstd_l0=std(Eccentricity_l0);
        ell_l0=mean(Ellipticity_l0);
        ellstd_l0=std(Ellipticity_l0);
        Ecc_l0(i+1,j)=ecc_l0;Eccstd_l0(i+1,j)=eccstd_l0;
        Ell_l0(i+1,j)=ell_l0;Ellstd_l0(i+1,j)=ellstd_l0;
        SI_l0(i+1,j)=(mean(Imax_l0.^2)-mean(Imax_l0)^2)/(mean(Imax_l0)^2);
    end
end
%% Plot
shape=["s","d","o"];
fs=19;
figure(1)
for i=1:3
    errorbar(L0_t,Ecc_L0(i,:),Eccstd_L0(i,:),shape(i),'MarkerSize',8,'LineWidth',1.5);
    hold on;
end 
hold off;
legend('m=0 modified','m=1 modified','m=2 modified','Location','Best');
xlabel('Outer Scale L_{0} (m)');
ylabel('Eccentricity');
axis([0 1 0.1 0.8]);
set(gca,'FontSize',fs);
figure(2)
for i=1:3
    errorbar(L0_t,Ell_L0(i,:),Ellstd_L0(i,:),shape(i),'MarkerSize',8,'LineWidth',1.5);
    hold on;
end 
hold off;
legend('m=0 modified','m=1 modified','m=2 modified','Location','Best');
xlabel('Outer Scale L_{0} (m)');
ylabel('Ellipticity');
axis([0 1 0.7 1]);
set(gca,'FontSize',fs);
figure(3)
color=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];
for i=1:3
    plot(L0_t,SI_L0(i,:),'o','MarkerSize',10,'MarkerFaceColor',color(i,:),'MarkerEdgeColor',color(i,:));hold on;
end
hold off;
legend('m=0 modified','m=1 modified','m=2 modified','Location','Best');
xlabel('Outer Scale L_{0} (m)');
ylabel('Scintillation Index');
axis([0 1 0 0.3]);
set(gca,'FontSize',fs);
figure(4)
for i=1:3
    errorbar(l0_t,Ecc_l0(i,:),Eccstd_l0(i,:),shape(i),'MarkerSize',8,'LineWidth',1.5);
    hold on;
end 
hold off;
legend('m=0 modified','m=1 modified','m=2 modified','Location','Best');
xlabel('Inner Scale l_{0} (m)');
ylabel('Eccentricity');
axis([0 0.01 0.1 0.8]);
set(gca,'FontSize',fs);
figure(5)
for i=1:3
    errorbar(l0_t,Ell_l0(i,:),Ellstd_l0(i,:),shape(i),'MarkerSize',8,'LineWidth',1.5);
    hold on;
end 
hold off;
legend('m=0 modified','m=1 modified','m=2 modified','Location','Best');
xlabel('Inner Scale l_{0} (m)');
ylabel('Ellipticity');
axis([0 0.01 0.7 1]);
set(gca,'FontSize',fs);
figure(6)
color=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];
for i=1:3
    plot(l0_t,SI_l0(i,:),'o','MarkerSize',10,'MarkerFaceColor',color(i,:),'MarkerEdgeColor',color(i,:));hold on;
end
hold off;
legend('m=0 modified','m=1 modified','m=2 modified','Location','Best');
xlabel('Inner Scale l_{0} (m)');
ylabel('Scintillation Index');
axis([0 0.01 0 0.3]);
set(gca,'FontSize',fs);
%%
save('simulation_result_all_1009_scale_length_2E-10.mat','Ecc_L0','Eccstd_L0','Ell_L0','Ellstd_L0','SI_L0','l0_t','L0_t',"Ecc_l0",'Ell_l0','SI_l0','Eccstd_l0',"Ellstd_l0");
