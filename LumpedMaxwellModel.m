clear all;

%% Set Parameters
tf = 24*200*3600;               %total integration time (in seconds)
nt = 24*1200;                   %number of time steps
dt = tf/nt;                     %time step length

w_tide1 = 24*3600*0.5;          %period of first tidal component
w_tide2 = 24*3600*0.518;        %period of secondst tidal component

length_scale = 20e3;            %Length scale between two stations to calculate strain

%Physical params
E = 9e9;                        %Elastic modulus
A_Glen=7e-25;                   %Nye-Glen rate factor
n=3;                            %Nye-Glen rate exponent

%choice params
alpha = 1.54;                   %Exponent of buttressing stress w.r.t tidal height
tidal_amp = 3;                  %Amplitude of combined tides
sigma_b0 = 75e3;                %Amplitude of tidal changes in buttressing stress
sigma_h0 = 55e3;                %Amplitude of tidal changes in hydrostatic stress
sigma_other = 160e3;            %Total background stresses

beta = sigma_b0/sigma_h0;       %Calculate beta non-dimensional parameter
gamma = sigma_other/sigma_h0;   %Calculate gamma non-dimensional parameter

%% Spin up integration
sigmas = [];
hs = [];
disps = [];
vels = [];
srs = [];

disp = 0;
ldt = 3600*24/4;                %Time step for spinup
for t = 0:ldt:ldt*5e3
    h_star = 0.5.*(sin(2*pi*t/w_tide1) + sin(2*pi*t/w_tide2));                      %normalized tidal height
    sigma = sigma_h0*(beta*(2*(0.5*(h_star+1)).^alpha - 1) - h_star);               %calculate stress

    dhs_dt = 0.5.*((2*pi/w_tide1).*cos(2*pi*t/w_tide1) + (2*pi/w_tide2).*cos(2*pi*t/w_tide2)); %calculate derivative of tidal height
    dsigma_dt = sigma_h0.*dhs_dt.*(alpha*beta*((0.5*(h_star+1)).^(alpha-1)) - 1);              %total stress derivative

    sigma_eff = sigma_h0.*sqrt(0.5*(beta.^2)*(2*(0.5*(h_star+1)).^(alpha) - 1)^2 + 0.5*(h_star.^2) + gamma^2);    %invariant of stress tensor

    velocity = length_scale*((1/(E)).*dsigma_dt + (2/3)*A_Glen.*(sigma_eff.^(n-1)).*sigma);     %calculate relative velocity between stations
    disp  = disp + velocity*dt;                                                                 %calculate displacement

end

%% Integration

for t = 0:dt:tf
    h_star = 0.5.*(sin(2*pi*t/w_tide1) + sin(2*pi*t/w_tide2));                      %normalized tidal height
    sigma = sigma_h0*(beta*(2*(0.5*(h_star+1)).^alpha - 1) - h_star);               %calculate stress

    dhs_dt = 0.5.*((2*pi/w_tide1).*cos(2*pi*t/w_tide1) + (2*pi/w_tide2).*cos(2*pi*t/w_tide2)); %calculate derivative of tidal height
    dsigma_dt = sigma_h0.*dhs_dt.*(alpha*beta*((0.5*(h_star+1)).^(alpha-1)) - 1);              %total stress derivative

    sigma_eff = sigma_h0.*sqrt(0.5*(beta.^2)*(2*(0.5*(h_star+1)).^(alpha) - 1)^2 + 0.5*(h_star.^2) + gamma^2);    %invariant of stress tensor

    velocity = length_scale*((1/(E)).*dsigma_dt + (2/3)*A_Glen.*(sigma_eff.^(n-1)).*sigma);     %calculate relative velocity between stations
    strain_rate = velocity./length_scale;
    disp  = disp + velocity*dt;                                                                 %calculate displacement
    
    %store outputs
    hs = [hs;h_star];
    sigmas = [sigmas;sigma];
    disps = [disps;disp];
    vels = [vels;velocity];
    srs = [srs;strain_rate];
end

ts=linspace(0,tf,nt+1)'/(24*3600);                      %time variable
disps_detrend = disps-polyval(polyfit(ts,disps,1),ts);  %detrend displacement

%% Calculate various characteristics of model out
T = 1/((ts(2)-ts(1))*24*60);
Fs = 1/T;
L = length(hs);
NFFT = 2^nextpow2(L);           
f = Fs/2*linspace(0,1,NFFT/2+1);
fft_Z = fft(hs,NFFT)/L;                 %Calculate fft of tidal height
fft_x20 = fft(disps_detrend,NFFT)/L;    %Calculate fft of detrended displacement
PSD = 2*abs(fft_Z(1:NFFT/2+1));
idx_primarytide = find(PSD==max(PSD));
phase_lag_Zx20 = 180*angle(fft_Z(find(PSD==max(PSD)))/fft_x20(find(PSD==max(PSD))))/pi; %phase lag between tidal height and displacement
    
%% Plot
figure;
[AX,H1,H2] = plotyy(ts-36,100*disps_detrend,ts-36,tidal_amp.*hs);hold on
xlabel('Time (days)','fontsize',22)
set(get(AX(2),'Ylabel'),'String','Tidal Height (m)','fontsize',22) 
ylabh = get(AX(2),'Ylabel');
set(get(AX(1),'Ylabel'),'String','Detrended Displacement (cm)','fontsize',22) 
set(AX(2),'xlim',[0 15],'ylim',[-3 3],'fontsize',22,'linewidth',2,'XTick',0:5:15,'Ytick',-3:1.5:3)
set(AX(1),'xlim',[0 15],'ylim',[-40 40],'fontsize',22,'linewidth',2,'XTick',0:5:15,'Ytick',-40:20:40)
set(H1,'linewidth',2);set(H2,'linewidth',2);set(gca,'fontsize',22)
