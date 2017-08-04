clear all
close all

%% Rutford
[daynr00,RotX00,RotY00,RotZ00]=ReadLocalXYZ('./Rutford/R+00.kin');
[daynr20,RotX20,RotY20,RotZ20]=ReadLocalXYZ('./Rutford/R-20.kin');

daynr20 = daynr20 - 731933;
daynr00 = daynr00 - 731933;
station_sep = 20e3;

% Interpolate to common times
time = (max([min(daynr00),min(daynr20)]):(5/(60*24)):min([max(daynr00),max(daynr20)]))';

x20_interp = interp1(daynr20,RotX20,time);
x00_interp = interp1(daynr00,RotX00,time);
z20_interp = interp1(daynr20,RotZ20,time);

% Make displacement

Lx_020 = x20_interp - x00_interp;

% Detrend

Px1=polyfit(time,x20_interp,1);
Detrendx20=x20_interp-polyval(Px1,time);

Px2=polyfit(time,x00_interp,1);
Detrendx00=x00_interp-polyval(Px2,time);

Px=polyfit(time,Lx_020,1);
DetrendLx=Lx_020-polyval(Px,time);

Pz=polyfit(time,z20_interp,1);
DetrendZ=z20_interp-polyval(Pz,time);

% Despike here
DS_thresh = 1.5;

x20_stdev = std(Detrendx20);
x00_stdev = std(Detrendx00);
Z_stdev = std(DetrendZ);

Despikex20 = Detrendx20;Despikex00 = Detrendx00;DespikeLx = DetrendLx;DespikeZ = DetrendZ;

DetrendLx_hf = DetrendLx - smooth(DetrendLx,300);
Lxhf_std = std(DetrendLx_hf);

Despikex20(abs(Detrendx20)>DS_thresh*x20_stdev | abs(Detrendx00)>DS_thresh*x00_stdev | abs(DetrendZ)>2*DS_thresh*Z_stdev | abs(DetrendLx_hf)>1*Lxhf_std) = [];
Despikex00(abs(Detrendx20)>DS_thresh*x20_stdev | abs(Detrendx00)>DS_thresh*x00_stdev | abs(DetrendZ)>2*DS_thresh*Z_stdev | abs(DetrendLx_hf)>1*Lxhf_std) = [];
time(abs(Detrendx20)>DS_thresh*x20_stdev | abs(Detrendx00)>DS_thresh*x00_stdev | abs(DetrendZ)>2*DS_thresh*Z_stdev | abs(DetrendLx_hf)>1*Lxhf_std) = [];
DespikeLx(abs(Detrendx20)>DS_thresh*x20_stdev | abs(Detrendx00)>DS_thresh*x00_stdev | abs(DetrendZ)>2*DS_thresh*Z_stdev | abs(DetrendLx_hf)>1*Lxhf_std) = [];
DespikeZ(abs(Detrendx20)>DS_thresh*x20_stdev | abs(Detrendx00)>DS_thresh*x00_stdev | abs(DetrendZ)>2*DS_thresh*Z_stdev | abs(DetrendLx_hf)>1*Lxhf_std) = [];


N = 6;

windowWidth = int16(N);
halfWidth = windowWidth / 2;

gaussFilter = gausswin(N);
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.

x20_dt_gf = conv(Despikex20, gaussFilter);
x00_dt_gf = conv(Despikex00, gaussFilter);
lx_dt_gf = conv(DespikeLx, gaussFilter);
Z_dt_gf = conv(DespikeZ, gaussFilter);

DespikeLx = lx_dt_gf(halfWidth:end-halfWidth);
DespikeZ = Z_dt_gf(halfWidth:end-halfWidth);

%Phase lags
T = 1/12;
Fs = 12;
L = length(DespikeZ);
NFFT = 2^nextpow2(L);
f = Fs/2*linspace(0,1,NFFT/2+1);
fft_Z = fft(DespikeZ,NFFT)/L;
fft_Lx = fft(DespikeLx,NFFT)/L;
fft_x20 = fft(Despikex20,NFFT)/L;
fft_x00 = fft(Despikex00,NFFT)/L;

PSD = 2*abs(fft_Z(1:NFFT/2+1));
idx_primarytide = find(PSD==max(PSD));

phase_lag_ZLx = 180*angle(fft_Z(find(PSD==max(PSD)))/fft_Lx(find(PSD==max(PSD))))/pi;
phase_lag_Zx20 = 180*angle(fft_Z(find(PSD==max(PSD)))/fft_x20(find(PSD==max(PSD))))/pi;
phase_lag_Zx00 = 180*angle(fft_Z(find(PSD==max(PSD)))/fft_x00(find(PSD==max(PSD))))/pi;


% DespikeLx= DespikeLx-smooth(DespikeLx,480);
% DespikeLx = Despikex20-smooth(Despikex20,480);
% Plots
figure(1);set(1,'units','pixels','position',[0 0 1800 1100])

subplot(3,2,1)
[AX,H1,H2] = plotyy(time-time(1)-38.5,100*Despikex20,time-time(1)-38.5,DespikeZ);hold on
xlabel('Time (days)','fontsize',22);
set(get(AX(2),'Ylabel'),'String','Tidal Height (m)','fontsize',22,'Rotation',270)
ylabh = get(AX(2),'Ylabel');
set(ylabh,'Position',get(ylabh,'Position') - [4.8 0 0])
set(get(AX(1),'Ylabel'),'String','Detrended Displacement (cm)','fontsize',22) 
set(AX(2),'fontsize',22,'linewidth',1,'Ylim',[-3 3],'Xlim',[0 15],'YTick',-3:1.5:3,'XTick',[0:5:15])
set(AX(1),'fontsize',22,'linewidth',1,'Ylim',[-40 40],'Xlim',[0 15],'YTick',-40:20:40,'XTick',[0:5:15])
set(H1,'linewidth',3);set(H2,'linewidth',2);set(gca,'fontsize',22)
text(0.01,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'Interpreter','LaTeX')
% text(0.05,0.12,['$M_2$ phase lag = ' num2str(round(phase_lag_Zx20*10)/10) '$^{\circ}$'],'Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',20,'Interpreter','LaTeX','BackgroundColor',[1 1 1],'Linewidth',1,'EdgeColor','k')
title('Observations: Rutford Ice Stream (R-20)','fontsize',22,'Interpreter','LaTeX')

tidal_height = DespikeZ-min(DespikeZ);
f = fittype('a+(b*(x^c))');
options = fitoptions('gauss2','Lower', [-0.2 -0.2 0],'Upper',[0.2 0.2 2],'TolFun',1e-9);
[fit1,gof,fitinfo] = fit(tidal_height,DespikeLx,f,options);
hs = linspace(min(tidal_height),max(tidal_height),100);
mid = coeffvalues(fit1);
cis = confint(fit1,0.95);

% [Pxx,F] = periodogram(Detrendx20,[],length(Detrendx20),3600*24/300);
% axes('Position',[.335 .735 .12 .05])
% box on
% plot((1./F),10*log10(Pxx),'k','linewidth',4);xlim([0 20]);ylim([-60 4])
% set(gca,'fontsize',16,'XTick',[0,20],'Ytick',[-60 -30 0]);ylabel('PSD','fontsize',16);%xlabel('Period (days)','fontsize',12);
% text(.25,-0.05,'Period (days)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',16)
% text(0.01,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'Interpreter','LaTeX')

subplot(3,2,3)
hf=plot(DespikeZ,(1e6/station_sep)*DespikeLx,'.','Color',[0.5 0.5 0.5],'markersize',5);hold on
set(get(get(hf,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
plot(hs+min(DespikeZ),(1e6/station_sep)*(mid(2).*(hs.^mid(3))+mid(1)),'k','linewidth',3)
plot(hs+min(DespikeZ),(1e6/station_sep)*(cis(1,2).*(hs.^cis(1,3))+cis(1,1)),'k--','linewidth',3)
plot(hs+min(DespikeZ),(1e6/station_sep)*(cis(2,2).*(hs.^cis(2,3))+cis(2,1)),'k--','linewidth',3)
xlabel('Tidal Height (m)','fontsize',22)
ylabel('Strain (10^{-6})','fontsize',22) 
axis([ -3 3 -5 7.5])
set(gca,'fontsize',24,'YTick',-5:5:5)
text(0.01,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'Interpreter','LaTeX')
% text(0.05,0.12,['$M_2$ phase lag = ' num2str(round(phase_lag_ZLx*10)/10) '$^{\circ}$'],'Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',20,'Interpreter','LaTeX','BackgroundColor',[1 1 1],'Linewidth',1,'EdgeColor','k')
title('Observations: Rutford Ice Stream (R-20 - R+00)','fontsize',22,'Interpreter','LaTeX')
hl=legend(['$\epsilon = ' num2str(round((1e6/station_sep)*mid(1),2)) ' + ' num2str(round((1e6/station_sep)*mid(2),2)) '(' num2str(round(abs(min(DespikeZ)),2)) '+h)^{' num2str(round(mid(3),2)) '}$'],'Location','SouthEast')
set(hl,'Interpreter','LaTeX','fontsize',20)


%% Bindschadler
site1 = 'DFLT';
site2 = 'D010';
station_sep = 20e3;

C_scale=[5;2.5;1.5];%conservative scales for e,n,u

file1 = ['./Bindschadler_data_fromdoug/smap.' site1 '.enu'];
file2 = ['./Bindschadler_data_fromdoug/smap.' site2 '.enu'];


load(file1); % variable is now called smap
load(file2); % variable is now called smap

%compute 2d displacement 
displ_DFLT=sqrt(smap_DFLT(1:end,2).^2+smap_DFLT(1:end,3).^2);
displ_D010=sqrt(smap_D010(1:end,2).^2+smap_D010(1:end,3).^2);

daynr00 = smap_D010(:,1);
daynr20 = smap_DFLT(:,1);

RotX00 = displ_D010;    %name mapping
RotX20 = displ_DFLT;    %name mapping

RotZ00 = smap_D010(:,4);
RotZ20 = smap_DFLT(:,4);

%interpolate
time = (max([min(daynr00),min(daynr20)]):(5/(60*24)):min([max(daynr00),max(daynr20)]))';

x20_interp = interp1(daynr20,RotX20,time);
x00_interp = interp1(daynr00,RotX00,time);
z20_interp = interp1(daynr20,RotZ20,time);

%make displacement
Lx_020 = x20_interp - x00_interp;

% Detrend
Px1=polyfit(time,x20_interp,1);
Detrendx20=x20_interp-polyval(Px1,time);

Px2=polyfit(time,x00_interp,1);
Detrendx00=x00_interp-polyval(Px2,time);

Px=polyfit(time,Lx_020,1);
DetrendLx=Lx_020-polyval(Px,time);

Pz=polyfit(time,z20_interp,1);
DetrendZ=z20_interp-polyval(Pz,time);

% Despike
DS_thresh = 2;

x20_stdev = std(Detrendx20);
x00_stdev = std(Detrendx00);
Z_stdev = std(DetrendZ);

Despikex20 = Detrendx20;Despikex00 = Detrendx00;DespikeLx = DetrendLx;DespikeZ = DetrendZ;

Despikex20(abs(Detrendx20)>DS_thresh*x20_stdev | abs(Detrendx00)>DS_thresh*x00_stdev | abs(DetrendZ)>2*DS_thresh*Z_stdev) = [];
Despikex00(abs(Detrendx20)>DS_thresh*x20_stdev | abs(Detrendx00)>DS_thresh*x00_stdev | abs(DetrendZ)>2*DS_thresh*Z_stdev) = [];
time(abs(Detrendx20)>DS_thresh*x20_stdev | abs(Detrendx00)>DS_thresh*x00_stdev | abs(DetrendZ)>2*DS_thresh*Z_stdev) = [];
DespikeLx(abs(Detrendx20)>DS_thresh*x20_stdev | abs(Detrendx00)>DS_thresh*x00_stdev | abs(DetrendZ)>2*DS_thresh*Z_stdev) = [];
DespikeZ(abs(Detrendx20)>DS_thresh*x20_stdev | abs(Detrendx00)>DS_thresh*x00_stdev | abs(DetrendZ)>2*DS_thresh*Z_stdev) = [];

N = 6;

windowWidth = int16(N);
halfWidth = windowWidth / 2;

gaussFilter = gausswin(N);
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.

x20_dt_gf = conv(Despikex20, gaussFilter);
x00_dt_gf = conv(Despikex00, gaussFilter);
lx_dt_gf = conv(DespikeLx, gaussFilter);
Z_dt_gf = conv(DespikeZ, gaussFilter);

DespikeLx = lx_dt_gf(halfWidth:end-halfWidth);
DespikeZ = Z_dt_gf(halfWidth:end-halfWidth);

%phase lags
T = 1/12;
Fs = 12;
L = length(DespikeZ);
NFFT = 2^nextpow2(L);
f = Fs/2*linspace(0,1,NFFT/2+1);
fft_Z = fft(DespikeZ,NFFT)/L;
fft_Lx = fft(DespikeLx,NFFT)/L;
fft_x20 = fft(Despikex20,NFFT)/L;
fft_x00 = fft(Despikex00,NFFT)/L;

PSD = 2*abs(fft_Z(1:NFFT/2+1));
idx_primarytide = find(PSD==max(PSD));

phase_lag_ZLx = 180*angle(fft_Z(find(PSD==max(PSD)))/fft_Lx(find(PSD==max(PSD))))/pi;
phase_lag_Zx20 = 180*angle(fft_Z(find(PSD==max(PSD)))/fft_x20(find(PSD==max(PSD))))/pi;
phase_lag_Zx00 = 180*angle(fft_Z(find(PSD==max(PSD)))/fft_x00(find(PSD==max(PSD))))/pi;

% DespikeLx= DespikeLx-smooth(DespikeLx,480);
% DespikeLx = Despikex20-smooth(Despikex20,240);

% Plots
figure(1);

subplot(3,2,2)
[AX,H1,H2] = plotyy(time-time(1)-29,100*Despikex20,time-time(1)-29,DespikeZ);hold on
xlabel('Time (days)','fontsize',22);
set(get(AX(2),'Ylabel'),'String','Tidal Height (m)','fontsize',22,'Rotation',270)
ylabh = get(AX(2),'Ylabel');
set(ylabh,'Position',get(ylabh,'Position') - [35.3 0 0]) 
set(get(AX(1),'Ylabel'),'String','Detrended Displacement (cm)','fontsize',22) 
set(AX(2),'fontsize',22,'linewidth',2,'Ylim',[-1.6 1.6],'Xlim',[0 15],'YTick',-1.6:0.8:1.6,'XTick',[0:5:15])
set(AX(1),'fontsize',22,'linewidth',2,'Ylim',[-12 12],'Xlim',[0 15],'YTick',-12:6:12,'XTick',[0:5:15])
set(H1,'linewidth',3);set(H2,'linewidth',2);set(gca,'fontsize',22)
text(0.01,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'Interpreter','LaTeX')
% text(0.05,0.125,['$O_1$ phase lag = ' num2str(360+round(phase_lag_Zx20*10)/10) '$^{\circ}$'],'Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',20,'Interpreter','LaTeX','BackgroundColor',[1 1 1],'Linewidth',1,'EdgeColor','k')
title('Observations: Bindschadler Ice Stream (DFLT)','fontsize',22,'Interpreter','LaTeX')

% flip_point = -1.2;
flip_point = min(-DespikeZ);
tidal_height = -DespikeZ-flip_point;
% DespikeLx(tidal_height<0)=[];DespikeZ(tidal_height<0)=[];tidal_height(tidal_height<0)=[];
f = fittype('a+(b*(x^c))');
options = fitoptions('gauss2','Lower', [-0.2 -0.2 0],'Upper',[0.2 0.2 5]);
[fit1,gof,fitinfo] = fit(tidal_height,DespikeLx,f,options);
hs = linspace(min(tidal_height),max(tidal_height),100);
mid = coeffvalues(fit1);
cis = confint(fit1,0.95);

% [Pxx,F] = periodogram(Detrendx20,[],length(DetrendLx),3600*24/300);
% % subplot(5,2,6)
% axes('Position',[.777 .735 .12 .05])
% box on
% plot((1./F),10*log10(Pxx),'k','linewidth',4);xlim([0 20]);ylim([-60 4])
% set(gca,'fontsize',16,'XTick',[0,20],'Ytick',[-60 -30 0]);ylabel('PSD','fontsize',16);%xlabel('Period (days)','fontsize',12);
% text(.25,-0.05,'Period (days)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',16)
% plot((1./F),10*log10(Pxx),'k','linewidth',4);xlim([0 20]);ylim([-60 2])
% set(gca,'fontsize',22);xlabel('Period (days)','fontsize',22);ylabel('log(Power)','fontsize',22)
% text(0.01,0.98,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'Interpreter','LaTeX')

subplot(3,2,4)
hf=plot(DespikeZ,1e6*DespikeLx./station_sep,'.','Color',[0.5 0.5 0.5],'markersize',5);hold on
set(get(get(hf,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
plot(hs+flip_point,1e6*(1/station_sep)*(mid(2).*((-2*flip_point-hs).^mid(3))+mid(1)),'k','linewidth',3)
plot(hs+flip_point,1e6*(1/station_sep)*(cis(1,2).*((-2*flip_point-hs).^cis(1,3))+cis(1,1)),'k--','linewidth',3)
plot(hs+flip_point,1e6*(1/station_sep)*(cis(2,2).*((-2*flip_point-hs).^cis(2,3))+cis(2,1)),'k--','linewidth',3)
xlabel('Tidal Height (m)','fontsize',22)
ylabel('Strain (10^{-6})','fontsize',22) 
axis([ -1.5 1.5 -0.05*1e6/station_sep 0.1*1e6/station_sep])
set(gca,'fontsize',24)
title('Observations: Bindschadler Ice Stream (DFLT - D010)','fontsize',22,'Interpreter','LaTeX')
text(0.01,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'Interpreter','LaTeX')

% text(0.05,0.12,['$O_1$ phase lag = ' num2str(360+round(phase_lag_ZLx*10)/10) '$^{\circ}$'],'Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',20,'Interpreter','LaTeX','BackgroundColor',[1 1 1],'Linewidth',1,'EdgeColor','k')
hl2=legend(['$\epsilon = ' num2str(round((1e6/station_sep)*mid(1),2)) ' + ' num2str((1e6/station_sep)*round(mid(2),3)) '(' num2str(round(-2*flip_point,2)) '-h)^{' num2str(round(mid(3),2)) '}$' ],'Location','NorthEast')
set(hl2,'Interpreter','LaTeX','fontsize',20)
% title(['\alpha = ' num2str(mid(3),3) ';2\sigma CI: ' num2str(cis(1,3),3) '<\alpha<' num2str(cis(2,3),3) '; r^2 =' num2str(gof.rsquare)],'fontsize',20)

%% Simple model - Rutford limit

clear all;

sigmas = [];
hs = [];
disps = [];
vels = [];
srs = [];

%Basic params
tf = 24*200*3600;
nt = 24*1200;
dt = tf/nt;

w_tide1 = 24*3600*0.5;
w_tide2 = 24*3600*0.518;

length_scale = 20e3;

%Physical params
E = 9e9;
A_Glen=7e-25;
n=3;

%choice params
alpha = 1.54;
tidal_amp = 3;
sigma_b0 = 75e3;
sigma_h0 = 55e3;
sigma_other = 160e3;

beta = sigma_b0/sigma_h0;
gamma = sigma_other/sigma_h0;

disp = 0;
ldt = 3600*24/4;
for t = 0:ldt:ldt*5e3
    h_star = 0.5.*(sin(2*pi*t/w_tide1) + sin(2*pi*t/w_tide2));  %normalized tidal height
    sigma = sigma_h0*(beta*(2*(0.5*(h_star+1)).^alpha - 1) - h_star);

    dhs_dt = 0.5.*((2*pi/w_tide1).*cos(2*pi*t/w_tide1) + (2*pi/w_tide2).*cos(2*pi*t/w_tide2)); 
    dsigma_dt = sigma_h0.*dhs_dt.*(alpha*beta*((0.5*(h_star+1)).^(alpha-1)) - 1); %total elastic stress derivative

    sigma_eff = sigma_h0.*sqrt(0.5*(beta.^2)*(2*(0.5*(h_star+1)).^(alpha) - 1)^2 + 0.5*(h_star.^2) + gamma^2);    %invariant of stress tensor

    velocity = length_scale*((1/(E)).*dsigma_dt + (2/3)*A_Glen.*(sigma_eff.^(n-1)).*sigma);
    disp  = disp + velocity*dt;

end


for t = 0:dt:tf
    h_star = 0.5.*(sin(2*pi*t/w_tide1) + sin(2*pi*t/w_tide2));  %normalized tidal height
    sigma = sigma_h0*(beta*(2*(0.5*(h_star+1)).^alpha - 1) - h_star);

    dhs_dt = 0.5.*((2*pi/w_tide1).*cos(2*pi*t/w_tide1) + (2*pi/w_tide2).*cos(2*pi*t/w_tide2)); 
    dsigma_dt = sigma_h0.*dhs_dt.*(alpha*beta*((0.5*(h_star+1)).^(alpha-1)) - 1); %total elastic stress derivative

    sigma_eff = sigma_h0.*sqrt(0.5*(beta.^2)*(2*(0.5*(h_star+1)).^(alpha) - 1)^2 + 0.5*(h_star.^2) + gamma^2);    %invariant of stress tensor

    velocity = length_scale*((1/(E)).*dsigma_dt + (2/3)*A_Glen.*(sigma_eff.^(n-1)).*sigma);
    strain_rate = velocity./length_scale;
    disp  = disp + velocity*dt;
    
    hs = [hs;h_star];
    sigmas = [sigmas;sigma];
    disps = [disps;disp];
    vels = [vels;velocity];
    srs = [srs;strain_rate];
end

ts=linspace(0,tf,nt+1)'/(24*3600);
disps_detrend = disps-polyval(polyfit(ts,disps,1),ts);

%phase lags
T = 1/((ts(2)-ts(1))*24*60);
Fs = 1/T;
L = length(hs);
NFFT = 2^nextpow2(L);
f = Fs/2*linspace(0,1,NFFT/2+1);
fft_Z = fft(hs,NFFT)/L;
fft_x20 = fft(disps_detrend,NFFT)/L;
PSD = 2*abs(fft_Z(1:NFFT/2+1));
idx_primarytide = find(PSD==max(PSD));
phase_lag_Zx20 = 180*angle(fft_Z(find(PSD==max(PSD)))/fft_x20(find(PSD==max(PSD))))/pi

% figure(2);set(2,'units','pixels','position',[0 0 1802 802])
subplot(3,2,5)
[AX,H1,H2] = plotyy(ts-36,100*disps_detrend,ts-36,tidal_amp.*hs);hold on
xlabel('Time (days)','fontsize',22)
set(get(AX(2),'Ylabel'),'String','Tidal Height (m)','fontsize',22,'Rotation',-90) 
ylabh = get(AX(2),'Ylabel');
set(ylabh,'Position',get(ylabh,'Position') - [189.8 0 0])
set(get(AX(1),'Ylabel'),'String','Detrended Displacement (cm)','fontsize',22) 
set(AX(2),'xlim',[0 15],'ylim',[-3 3],'fontsize',22,'linewidth',2,'XTick',0:5:15,'Ytick',-3:1.5:3)
set(AX(1),'xlim',[0 15],'ylim',[-40 40],'fontsize',22,'linewidth',2,'XTick',0:5:15,'Ytick',-40:20:40)
set(H1,'linewidth',2);set(H2,'linewidth',2);set(gca,'fontsize',22)
text(0.01,0.98,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'Interpreter','LaTeX')
% text(0.05,0.125,['$M_2$ phase lag = ' num2str(round(phase_lag_Zx20*10)/10) '$^{\circ}$'],'Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',20,'Interpreter','LaTeX','BackgroundColor',[1 1 1],'Linewidth',1,'EdgeColor','k')
title(['Simulations: Rutford-like ($\alpha = ' num2str(alpha) '$; $\beta =' num2str(round(beta*10)/10) '$; $\gamma =' num2str(round(gamma*10)/10) '$)'],'Fontsize',22,'Interpreter','LaTeX')

% subplot(5,2,[8,10])
% plot(tidal_amp.*hs,100*disps_detrend,'.','Color',[0.5 0.5 0.5],'markersize',2)
% xlabel('Tidal Height (m)','fontsize',26)
% ylabel('Displacement (cm)','fontsize',26) 
% set(gca,'fontsize',24)

% [Pxx,F] = periodogram(disps_detrend./max(disps_detrend),[],length(disps_detrend),3600*24/dt);
% % subplot(3,2,5)
% axes('Position',[.335 .135 .12 .05])
% box on
% plot((1./F),10*log10(Pxx),'k','linewidth',4);xlim([0 20]);ylim([-70 30])
% set(gca,'fontsize',16,'XTick',[0,20],'Ytick',[-60 -30 0]);ylabel('PSD','fontsize',16);%xlabel('Period (days)','fontsize',12);
% text(.25,-0.05,'Period (days)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',16)

drawnow

%% Simple model - Bindschadler limit

% clear all;close all

sigmas = [];
hs = [];
disps = [];
vels = [];

%Basic params
tf = 24*200*3600;
nt = 24*1200;
dt = tf/nt;

w_tide1 = 24*3600*.9971;
w_tide2 = 24*3600*1.075;

length_scale = 25e3;

%Physical params
E = 9e9;
A_Glen=7e-25;
n=3;

%choice params
alpha = 2.53;
tidal_amp = 1.6;
sigma_b0 = 10e3;
sigma_h0 = 30e3;
sigma_other = 120e3;

beta = sigma_b0/sigma_h0;
gamma = sigma_other/sigma_h0;

%integrate to SS
disp = 0;
ldt = 3600*24/4;
for t = 0:ldt:ldt*5e3
    h_star = 0.5.*(sin(2*pi*t/w_tide1) + sin(2*pi*t/w_tide2));  %normalized tidal height
    sigma = sigma_h0*(beta*(2*(0.5*(h_star+1)).^alpha - 1) - h_star);

    dhs_dt = 0.5.*((2*pi/w_tide1).*cos(2*pi*t/w_tide1) + (2*pi/w_tide2).*cos(2*pi*t/w_tide2)); 
    dsigma_dt = sigma_h0.*dhs_dt.*(alpha*beta*((0.5*(h_star+1)).^(alpha-1)) - 1); %total elastic stress derivative

    sigma_eff = sigma_h0.*sqrt(0.5*(beta.^2)*(2*(0.5*(h_star+1)).^(alpha) - 1)^2 + 0.5*(h_star.^2) + gamma^2);    %invariant of stress tensor

    velocity = length_scale*((1/(E)).*dsigma_dt + (2/3)*A_Glen.*(sigma_eff.^(n-1)).*sigma);
    disp  = disp + velocity*dt;

end

%integrate to desired length
for t = 0:dt:tf
    h_star = 0.5.*(sin(2*pi*t/w_tide1) + sin(2*pi*t/w_tide2));  %normalized tidal height
    sigma = sigma_h0*(beta*(2*(0.5*(h_star+1)).^alpha - 1) - h_star);

    dhs_dt = 0.5.*((2*pi/w_tide1).*cos(2*pi*t/w_tide1) + (2*pi/w_tide2).*cos(2*pi*t/w_tide2)); 
    dsigma_dt = sigma_h0.*dhs_dt.*(alpha*beta*((0.5*(h_star+1)).^(alpha-1)) - 1); %total elastic stress derivative

    sigma_eff = sigma_h0.*sqrt(0.5*(beta.^2)*(2*(0.5*(h_star+1)).^(alpha) - 1)^2 + 0.5*(h_star.^2) + gamma^2);    %invariant of stress tensor

    velocity = length_scale*((1/(E)).*dsigma_dt + (2/3)*A_Glen.*(sigma_eff.^(n-1)).*sigma);
    disp  = disp + velocity*dt;
    
    hs = [hs;h_star];
    sigmas = [sigmas;sigma];
    disps = [disps;disp];
    vels = [vels;velocity];
end

ts=linspace(0,tf,nt+1)'/(24*3600);
disps_detrend = disps-polyval(polyfit(ts,disps,1),ts);

%phase lags
T = 1/((ts(2)-ts(1))*24*60);
Fs = 1/T;
L = length(hs);
NFFT = 2^nextpow2(L);
f = Fs/2*linspace(0,1,NFFT/2+1);
fft_Z = fft(hs,NFFT)/L;
fft_x20 = fft(disps_detrend,NFFT)/L;
PSD = 2*abs(fft_Z(1:NFFT/2+1));
idx_primarytide = find(PSD==max(PSD));
phase_lag_Zx20 = 180*angle(fft_Z(find(PSD==max(PSD)))/fft_x20(find(PSD==max(PSD))))/pi;

% figure(2);set(2,'units','pixels','position',[0 0 1802 802])
subplot(3,2,6)
[AX,H1,H2] = plotyy(ts-33,100*disps_detrend,ts-33,tidal_amp.*hs);hold on
xlabel('Time (days)','fontsize',22)
set(get(AX(2),'Ylabel'),'String','Tidal Height (m)','fontsize',22,'Rotation',-90) 
ylabh = get(AX(2),'Ylabel');
set(ylabh,'Position',get(ylabh,'Position') - [189.8 0 0]) 
set(get(AX(1),'Ylabel'),'String','Detrended Displacement (cm)','fontsize',22) 
set(AX(2),'xlim',[0 15],'fontsize',22,'linewidth',2,'XTick',0:5:15,'ylim',[-1.6 1.6],'YTick',[-1.6:0.8:1.6])
set(AX(1),'xlim',[0 15],'ylim',[-12 12],'fontsize',22,'linewidth',2,'XTick',0:5:15,'YTick',[-12:6:12])
set(H1,'linewidth',2);set(H2,'linewidth',2);set(gca,'fontsize',22)
text(0.01,0.98,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'Interpreter','LaTeX')
% text(0.05,0.125,['$O_1$ phase lag = ' num2str(360+round(phase_lag_Zx20*10)/10) '$^{\circ}$'],'Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',20,'Interpreter','LaTeX','BackgroundColor',[1 1 1],'Linewidth',1,'EdgeColor','k')
title(['Simulations: Bindschadler-like ($\alpha = ' num2str(alpha) '$; $\beta =' num2str(round(beta*10)/10) '$; $\gamma =' num2str(round(gamma*10)/10) '$)'],'Fontsize',22,'Interpreter','LaTeX')

% title(['Detrended Dispacement at GL; $\alpha = ' num2str(alpha) '$; $\beta_{max} =' num2str(force0^2/(force_amp^2)) '$'],'Fontsize',20,'Interpreter','LaTeX')

% [Pxx,F] = periodogram(disps_detrend./max(disps_detrend),[],length(disps_detrend),3600*24/dt);
% % subplot(3,2,6)
% axes('Position',[.773 .135 .12 .05])
% box on
% plot((1./F),10*log10(Pxx),'k','linewidth',4);xlim([0 20]);ylim([-70 30])
% set(gca,'fontsize',16,'XTick',[0,20],'Ytick',[-60 -30 0]);ylabel('PSD','fontsize',16);%xlabel('Period (days)','fontsize',12);
% text(.25,-0.05,'Period (days)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',16)
drawnow

% figure(10)
% plot(tidal_amp.*hs,100*disps_detrend,'.','Color',[0.5 0.5 0.5],'markersize',2)
% xlabel('Tidal Height (m)','fontsize',26)
% ylabel('Displacement (cm)','fontsize',26) 
% set(gca,'fontsize',24)
