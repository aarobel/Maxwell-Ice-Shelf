% %% Vary beta and gamma and look at relative PSD and phase lag
% 
clear all;close all

%Basic params
tf = 24*200*3600;
nt = 24*600;
dt = tf/nt;

w_tide1 = 24*3600*0.5;
w_tide2 = 24*3600*0.518;

length_scale = 20e3;

%Physical params
E = 9e9;
A_Glen=7e-25;
n=3;

%choice params
tidal_amp = 1;
sigma_h0 = 50e3;

% alphas = linspace(0.5,3,10);
% hf1 = 1-linspace(1^(1/4),0,20).^4;
% hf2 = 1+linspace(0,2^(1/4),20).^4;
% betas = [hf1,hf2(2:end)];
betas = linspace(0.01,2,40);
% beta = 2;
gammas = linspace(0,10,20);%sigma_other/sigma_h0;

q=-1;
for alpha = [1.5 2.5]
q=q+1;
%     for i = 1:length(gammas)
%         gamma=gammas(i);
%         for j = 1:length(betas)
%         beta = betas(j);
% 
% 
%         sigmas = [];
%         hs = [];
%         disps = [];
%         vels = [];
% 
%         %integrate to SS
%         disp = 0;
%         ldt = 3600*24/4;
%         for t = 0:ldt:ldt*5e3
%             h_star = 0.5.*(sin(2*pi*t/w_tide1) + sin(2*pi*t/w_tide2));  %normalized tidal height
%             sigma = sigma_h0*(beta*(2*(0.5*(h_star+1)).^alpha - 1) - h_star);
% 
%             dhs_dt = 0.5.*((2*pi/w_tide1).*cos(2*pi*t/w_tide1) + (2*pi/w_tide2).*cos(2*pi*t/w_tide2)); 
%             dsigma_dt = sigma_h0.*dhs_dt.*(alpha*beta*((0.5*(h_star+1)).^(alpha-1)) - 1); %total elastic stress derivative
% 
%                 sigma_eff = sigma_h0.*sqrt(0.5*(beta.^2)*(2*(0.5*(h_star+1)).^(alpha) - 1)^2 + 0.5*(h_star.^2) + gamma^2);    %invariant of stress tensor
% 
%             velocity = length_scale*((1/(E)).*dsigma_dt + (2/3)*A_Glen.*(sigma_eff.^(n-1)).*sigma);
%             disp  = disp + velocity*dt;
% 
%         end
% 
%         %integrate to desired length
%         for t = 0:dt:tf
%             h_star = 0.5.*(sin(2*pi*t/w_tide1) + sin(2*pi*t/w_tide2));  %normalized tidal height
%             sigma = sigma_h0*(beta*(2*(0.5*(h_star+1)).^alpha - 1) - h_star);
% 
%             dhs_dt = 0.5.*((2*pi/w_tide1).*cos(2*pi*t/w_tide1) + (2*pi/w_tide2).*cos(2*pi*t/w_tide2)); 
%             dsigma_dt = sigma_h0.*dhs_dt.*(alpha*beta*((0.5*(h_star+1)).^(alpha-1)) - 1); %total elastic stress derivative
% 
%                 sigma_eff = sigma_h0.*sqrt(0.5*(beta.^2)*(2*(0.5*(h_star+1)).^(alpha) - 1)^2 + 0.5*(h_star.^2) + gamma^2);    %invariant of stress tensor
% 
%             velocity = length_scale*((1/(E)).*dsigma_dt + (2/3)*A_Glen.*(sigma_eff.^(n-1)).*sigma);
%             disp  = disp + velocity*dt;
% 
%             hs = [hs;h_star];
%             sigmas = [sigmas;sigma];
%             disps = [disps;disp];
%             vels = [vels;velocity];
%         end
% 
%         ts=linspace(0,tf,nt+1)'/(24*3600);
%         disps_detrend = disps-polyval(polyfit(ts,disps,1),ts);
% 
%         %calculate relative pSD
%         [Pxx,F] = periodogram(disps_detrend./max(disps_detrend),[],length(disps_detrend),3600*24/dt);
%         PSD_relative(i,j) = Pxx(15)./Pxx(387);
% 
%         %calculate phase lag between tide and displacement
%         T = 1/(3600*24/dt);
%         Fs = 3600*24/dt;
%         L = length(hs);
%         NFFT = 2^nextpow2(L);
%         f = Fs/2*linspace(0,1,NFFT/2+1);
% 
%         fft_Z = fft(hs,NFFT)/L;
%         fft_Lx = fft(disps_detrend,NFFT)/L;
%         PSD = 2*abs(fft_Z(1:NFFT/2+1));
% 
%         idx_primarytide = find(PSD==max(PSD));
%         phase = angle(fft_Z(find(PSD==max(PSD)))/fft_Lx(find(PSD==max(PSD))));
% 
%         if(180*phase./pi<-30)
%             PHASE(i,j) = 2*pi+phase;
%         else
%             PHASE(i,j) = phase;
%         end 
%         BETA(i,j) = beta;
%         GAMMA(i,j) = gamma;
%         [i,j]
%         end
%     end

    if(q==0)
        load('ParamSpace_alpha1-5_v2.mat')
    else
        load('ParamSpace_alpha2-5_v2.mat')
    end

    figure(1);set(1,'units','pixels','position',[0 0 1102 802]);set(1,'Renderer','Painters')
    ax1 = subplot(2,2,q+1);
    map = brewermap(32,'RdBu');colormap(ax1,map)
    pcolor(BETA,GAMMA,180*PHASE./pi);hold on;
    shading('flat');h=colorbar;set(h,'fontsize',20);ylabel(h,'Phase ($^{\circ}$)','Interpreter','LaTeX');
    caxis([0 240])
%     [C,h] = contour(BETA,GAMMA,180*PHASE./pi,[30 30]);
%     h.LineWidth = 2;h.LineColor = [1 1 1];
    xlabel('$\beta$','fontsize',26,'Interpreter','LaTeX')
    ylabel('$\gamma$','fontsize',26,'Interpreter','LaTeX')
    set(gca,'fontsize',24)
    if(q==0);
%         save('ParamSpace_alpha1-5.mat')
        load('ParamSpace_alpha1-5_v2.mat')
        h=plot(1.4,2.9,'o','markersize',20);set(h,'MarkerEdgeColor','none','MarkerFaceColor','k')
        text(0.01,0.11,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'color','w','Interpreter','LaTeX')
    else
%         save('ParamSpace_alpha2-5.mat')
    load('ParamSpace_alpha2-5_v2.mat')
        text(0.01,0.11,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'color','w','Interpreter','LaTeX')
%         h=plot(0.3,4,'o','markersize',20);set(h,'MarkerEdgeColor','none','MarkerFaceColor','r')
    end
    title(['$\alpha = ' num2str(round(alpha*10)/10) '$'],'fontsize',25,'Interpreter','LaTeX')
    drawnow
    
    figure(1);%set(2,'units','pixels','position',[0 0 1102 402]);
    ax2=subplot(2,2,q+3);
    map = brewermap(32,'PuOr');colormap(ax2,map)
    pcolor(BETA,GAMMA,log10(PSD_relative));hold on;
    shading('flat');h=colorbar;set(h,'fontsize',20);ylabel(h,'$log(S_{M_{sf}}) - log(S_{M_2})$','Interpreter','LaTeX');caxis([-3 3])
%     [C,h] = contour(BETA,GAMMA,log10(PSD_relative),[1 1]);
%     h.LineWidth = 2;h.LineColor = [1 1 1];
    xlabel('$\beta$','fontsize',26,'Interpreter','LaTeX')
    ylabel('$\gamma$','fontsize',26,'Interpreter','LaTeX')
    set(gca,'fontsize',24)
    if(q==0);
        text(0.01,0.11,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'color','w','Interpreter','LaTeX')
        h=plot(1.4,2.9,'o','markersize',20);set(h,'MarkerEdgeColor','none','MarkerFaceColor','k')
    else
        text(0.01,0.11,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'color','w','Interpreter','LaTeX')
%         h=plot(0.3,4,'o','markersize',20);set(h,'MarkerEdgeColor','none','MarkerFaceColor','r')
    end
    title(['$\alpha = ' num2str(round(alpha*10)/10) '$'],'fontsize',25,'Interpreter','LaTeX')
    drawnow
    
end



% clear all;
% 
% %Basic params
% tf = 24*200*3600;
% nt = 24*600;
% dt = tf/nt;
% 
% w_tide1 = 24*3600*0.5;
% w_tide2 = 24*3600*0.518;
% 
% length_scale = 20e3;
% 
% %Physical params
% G = 2e9;
% A_Glen=1e-24;
% n=3;
% 
% %choice params
% tidal_amp = 1;
% sigma_h0 = 25e3;
% 
% alphas = linspace(0.5,3,10);
% % hf1 = 1-linspace(1^(1/3),0,15).^3;
% % hf2 = 1+linspace(0,2^(1/3),15).^3;
% % betas = [hf1,hf2(2:end)];
% % beta = 2;
% gammas = linspace(0,6,10);%sigma_other/sigma_h0;
% 
% q=-1;
% for beta = [0.5 2]
% q=q+1;
%     for i = 1:length(gammas)
%         gamma=gammas(i);
%         for j = 1:length(alphas)
%         alpha=alphas(j);
% 
% 
%         sigmas = [];
%         hs = [];
%         disps = [];
%         vels = [];
% 
%         %integrate to SS
%         disp = 0;
%         ldt = 3600*24/4;
%         for t = 0:ldt:ldt*5e3
%             h_star = 0.5.*(sin(2*pi*t/w_tide1) + sin(2*pi*t/w_tide2));  %normalized tidal height
%             sigma = sigma_h0*(beta*(2*(0.5*(h_star+1)).^alpha - 1) - h_star);
% 
%             dhs_dt = 0.5.*((2*pi/w_tide1).*cos(2*pi*t/w_tide1) + (2*pi/w_tide2).*cos(2*pi*t/w_tide2)); 
%             dsigma_dt = sigma_h0.*dhs_dt.*(alpha*beta*((0.5*(h_star+1)).^(alpha-1)) - 1); %total elastic stress derivative
% 
%             sigma_eff = sigma_h0.*sqrt((beta.^2)*(2*(0.5*(h_star+1)).^(alpha) - 1)^2 + (h_star.^2) + gamma^2);    %invariant of stress tensor
% 
%             velocity = length_scale*((1/(2*G)).*dsigma_dt + A_Glen.*(sigma_eff.^(n-1)).*sigma);
%             disp  = disp + velocity*dt;
% 
%         end
% 
%         %integrate to desired length
%         for t = 0:dt:tf
%             h_star = 0.5.*(sin(2*pi*t/w_tide1) + sin(2*pi*t/w_tide2));  %normalized tidal height
%             sigma = sigma_h0*(beta*(2*(0.5*(h_star+1)).^alpha - 1) - h_star);
% 
%             dhs_dt = 0.5.*((2*pi/w_tide1).*cos(2*pi*t/w_tide1) + (2*pi/w_tide2).*cos(2*pi*t/w_tide2)); 
%             dsigma_dt = sigma_h0.*dhs_dt.*(alpha*beta*((0.5*(h_star+1)).^(alpha-1)) - 1); %total elastic stress derivative
% 
%             sigma_eff = sigma_h0.*sqrt((beta.^2)*(2*(0.5*(h_star+1)).^(alpha) - 1)^2 + (h_star.^2) + gamma^2);    %invariant of stress tensor
% 
%             velocity = length_scale*((1/(2*G)).*dsigma_dt + A_Glen.*(sigma_eff.^(n-1)).*sigma);
%             disp  = disp + velocity*dt;
% 
%             hs = [hs;h_star];
%             sigmas = [sigmas;sigma];
%             disps = [disps;disp];
%             vels = [vels;velocity];
%         end
% 
%         ts=linspace(0,tf,nt+1)'/(24*3600);
%         disps_detrend = disps-polyval(polyfit(ts,disps,1),ts);
% 
%         %calculate relative pSD
%         [Pxx,F] = periodogram(disps_detrend./max(disps_detrend),[],length(disps_detrend),3600*24/dt);
%         PSD_relative(i,j) = Pxx(15)./Pxx(387);
% 
%         %calculate phase lag between tide and displacement
%         T = 1/(3600*24/dt);
%         Fs = 3600*24/dt;
%         L = length(hs);
%         NFFT = 2^nextpow2(L);
%         f = Fs/2*linspace(0,1,NFFT/2+1);
% 
%         fft_Z = fft(hs,NFFT)/L;
%         fft_Lx = fft(disps_detrend,NFFT)/L;
%         PSD = 2*abs(fft_Z(1:NFFT/2+1));
% 
%         idx_primarytide = find(PSD==max(PSD));
%         phase = angle(fft_Z(find(PSD==max(PSD)))/fft_Lx(find(PSD==max(PSD))));
% 
%         if(180*phase./pi<-30)
%             PHASE(i,j) = 2*pi+phase;
%         else
%             PHASE(i,j) = phase;
%         end 
%         ALPHA(i,j) = alpha;
%         GAMMA(i,j) = gamma;
%         [i,j]
%         end
%     end
% 
% %     figure(1);
% %     subplot(2,2,q+1)
% %     pcolor(BETA,GAMMA,PSD_relative);hold on
% %     shading('flat');h=colorbar;set(h,'fontsize',20);ylabel(h,'PSD Ratio');caxis([0 100])
% %     [C,h] = contour(BETA,GAMMA,PSD_relative,[1 1]);
% %     h.LineWidth = 2;h.LineColor = [1 1 1];
% %     xlabel('$\beta$','fontsize',26,'Interpreter','LaTeX')
% %     ylabel('$\gamma$','fontsize',26,'Interpreter','LaTeX')
% %     set(gca,'fontsize',24)
% %     title(['$\alpha = ' num2str(round(alpha*10)/10) '$'],'fontsize',26,'Interpreter','LaTeX')
%     
%     figure(1);set(1,'units','pixels','position',[0 0 1102 752])
%     subplot(2,2,q+1)
%     pcolor(ALPHA,GAMMA,180*PHASE./pi);hold on;
%     shading('flat');h=colorbar;set(h,'fontsize',20);ylabel(h,'Phase ($^{\circ}$)','Interpreter','LaTeX');
%     if(beta>1);caxis([0 40]);else;caxis([180 220]);end
% %     [C,h] = contour(BETA,GAMMA,180*PHASE./pi,[30 30]);
% %     h.LineWidth = 2;h.LineColor = [1 1 1];
%     xlabel('$\alpha$','fontsize',26,'Interpreter','LaTeX')
%     ylabel('$\gamma$','fontsize',26,'Interpreter','LaTeX')
%     set(gca,'fontsize',24)
%     if(q==0);
%         text(0.01,0.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'color','w','Interpreter','LaTeX')
%     else
%         text(0.01,0.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',25,'color','w','Interpreter','LaTeX')
%     end
%     title(['$\beta = ' num2str(beta) '$'],'fontsize',26,'Interpreter','LaTeX')
% end


