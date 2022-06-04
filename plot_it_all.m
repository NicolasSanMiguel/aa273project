% Nicolas San Miguel - May 2022
% AA 273 - Spring 2022


%% plot EKF
if plot_lonlat

    % load('pvtTRANSPLANT.mat')
    trange = 1:length(londata);
    figure;
    scatter(x_ekf(2,trange),x_ekf(1,trange),'.','LineWidth',2);
    hold on;
    % scatter(UBXpvt.lon,UBXpvt.lat,'*','LineWidth',2); % transplant data
    scatter(londata(trange),latdata(trange),'.','LineWidth',2);
    hold on;
    % for i = 1:length(x_ekf)
    % scatter(x_ekf(2,i),x_ekf(1,i),'.','LineWidth',5);
    % end

    xlabel('Longitude'); ylabel('Latitude');
    title('state x-y plane')
    legend('Estimated Position','True Position Meas.','other')
end

if plot_lonlat_ss
    %% plot EKF steady state
    % load('pvtTRANSPLANT.mat') % for validating true position while
    % debugging
    trange = 3600:3600*2;
    figure;
    scatter(x_ekf(2,trange),x_ekf(1,trange),'.','LineWidth',2);
    hold on;
    % scatter(UBXpvt.lon,UBXpvt.lat,'*','LineWidth',2); % transplant data
    scatter(londata(trange),latdata(trange),'.','LineWidth',2);
    hold on;
    xlabel('Longitude'); ylabel('Latitude');
    title('state x-y plane')
    legend('Estimated Position','True Position Meas.','other')
end
%% geoplotting
if plot_lonlat_geo
    figure
    trange = 3600:3600*2;
    geoscatter(latdata(trange),londata(trange),'*')
    hold on
    geoscatter(x_ekf(1,trange),x_ekf(2,trange),'o','filled')
    hold on
    legend('True','Estimated')
    geobasemap streets
end

if plot_states
    %% plot each state
    disp('here')
    figure
    subplot(3,1,1)
    %     trange = 11*3600:12*3600;
    plot(dateTvec(trange),londata(trange),'LineWidth',2)
    hold on
    plot(dateTvec(trange),x_ekf(2,trange),'LineWidth',1)
    if include_covariance
        perf = 1.96*sqrt(reshape(sig_ekf(2,2,trange),1,[]));
        patch([dateTvec(trange) fliplr(dateTvec(trange))], [x_ekf(2,trange)+perf fliplr(x_ekf(2,trange)-perf)],'k','FaceColor','#D95319','FaceAlpha',.2)
    end
    xlim([dateTvec(trange(1)) dateTvec(trange(end))])
    xlabel('Time'); ylabel('Longitude'); legend('True','Estimate','Covariance')


    subplot(3,1,2)
    plot(dateTvec(trange),latdata(trange),'LineWidth',2)
    hold on
    plot(dateTvec(trange),x_ekf(1,trange),'LineWidth',1)
    if include_covariance
        perf = 1.96*sqrt(reshape(sig_ekf(1,1,trange),1,[]));
        patch([dateTvec(trange) fliplr(dateTvec(trange))], [x_ekf(1,trange)+perf fliplr(x_ekf(1,trange)-perf)],'k','FaceColor','#D95319','FaceAlpha',.2)
    end
    xlim([dateTvec(trange(1)) dateTvec(trange(end))])
    xlabel('Time'); ylabel('Longitude'); legend('True','Estimate','Covariance')


    subplot(3,1,3)
    plot(dateTvec(trange),heightData(trange),'LineWidth',2)
    hold on
    plot(dateTvec(trange),x_ekf(3,trange),'LineWidth',2)
    if include_covariance
        perf = 1.96*sqrt(reshape(sig_ekf(3,3,trange),1,[]));
        patch([dateTvec(trange) fliplr(dateTvec(trange))], [x_ekf(3,trange)+perf fliplr(x_ekf(3,trange)-perf)],'k','FaceColor','#D95319','FaceAlpha',.2)
    end
    xlim([dateTvec(trange(1)) dateTvec(trange(end))])
    xlabel('Time'); ylabel('Height [m]'); legend('True','Estimate','Covariance')
end

if plot_states4
    %% plot each state
    disp('here')
    figure
    subplot(2,2,1)
    %     trange = 11*3600:12*3600;
    plot(dateTvec(trange),M_L1(trange),'LineWidth',2)
    hold on
    plot(dateTvec(trange),y(1,trange),'LineWidth',1)
    xlim([dateTvec(trange(1)) dateTvec(trange(end))])
    xlabel('Time'); ylabel('CN0_{Max} L1 [dB Hz]'); legend('True','Estimate','Covariance')

    subplot(2,2,2)
    %     trange = 11*3600:12*3600;
    plot(dateTvec(trange),M_L2(trange),'LineWidth',2)
    hold on
    plot(dateTvec(trange),y(2,trange),'LineWidth',1)
    xlim([dateTvec(trange(1)) dateTvec(trange(end))])
    xlabel('Time'); ylabel('CN0_{Max} L2 [dB Hz]'); legend('True','Estimate','Covariance')

    subplot(2,2,3)
    plot(dateTvec(trange),agcCntdata(1,trange),'LineWidth',2)
    hold on
    plot(dateTvec(trange),y(3,trange),'LineWidth',1)
    xlim([dateTvec(trange(1)) dateTvec(trange(end))])
    xlabel('Time'); ylabel('AGC L1'); legend('True','Estimate','Covariance')

    subplot(2,2,4)
    plot(dateTvec(trange),agcCntdata(2,trange),'LineWidth',2)
    hold on
    plot(dateTvec(trange),y(4,trange),'LineWidth',1)
    xlim([dateTvec(trange(1)) dateTvec(trange(end))])
    xlabel('Time'); ylabel('AGC L2'); legend('True','Estimate','Covariance')

%     subplot(3,2,5)
%     plot(dateTvec(trange),prResdata(2,trange),'LineWidth',2)
%     hold on
%     plot(dateTvec(trange),x_ekf(8,trange),'LineWidth',1)
%     if include_covariance
%         perf = 1.96*sqrt(reshape(sig_ekf(8,8,trange),1,[]));
%         patch([dateTvec(trange) fliplr(dateTvec(trange))], [x_ekf(8,trange)+perf fliplr(x_ekf(8,trange)-perf)],'k','FaceColor','#D95319','FaceAlpha',.2)
%     end
%     xlim([dateTvec(trange(1)) dateTvec(trange(end))])
%     xlabel('Time'); ylabel('prRes'); legend('True','Estimate','Covariance')

end


%% plot innovation vector
if plot_innov
    trange = 1:length(dateTvec);
    figure
    subplot(4,1,1)
    plot(dateTvec(trange),mrsInno(4,trange),'LineWidth',2)
    hold on
    % xlim([dateTvec(11*3600) dateTvec(12*3600)])
    xlabel('Time'); ylabel('CN0_{max} L1');
    title('Innovation Vector CN0max L1')

    subplot(4,1,2)
    plot(dateTvec(trange),mrsInno(5,trange),'LineWidth',2)
    hold on
    % xlim([dateTvec(11*3600) dateTvec(12*3600)])
    xlabel('Time'); ylabel('CN0_{max} L2');
    title('Innovation Vector CN0max L2')

    subplot(4,1,3)
    plot(dateTvec(trange),mrsInno(6,trange),'LineWidth',2)
    hold on
    % xlim([dateTvec(11*3600) dateTvec(12*3600)])
    xlabel('Time'); ylabel('AGC L1');
    title('Innovation Vector AGC L1')
    ylim([0,inf])

    subplot(4,1,4)
    plot(dateTvec(trange),mrsInno(7,trange),'LineWidth',2)
    hold on
    % xlim([dateTvec(11*3600) dateTvec(12*3600)])
    ylim([0,inf])
    xlabel('Time'); ylabel('AGC L2');
    title('Innovation Vector AGC L2')
end


%% helper functions
function errorEllipse(miu, sigma, color)
% this function takes mean and covariance matrix as parameters and plot
p = 0.95;     %probability
ellipse_const = -2 * log(1 - p);
r = sqrt(ellipse_const);
theta = 0:0.01:2 * pi;
w1 = r * cos(theta);
w2 = r * sin(theta);
w = [w1; w2];   % normalized coordinates
x = sqrtm(sigma) * w + miu * ones(1, size(w, 2));
plot(x(1,:), x(2,:), color);
end
