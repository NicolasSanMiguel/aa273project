% Nicolas San Miguel - May 2022
% AA 273 - Spring 2022

%% Raw Chi-squared
% mrsInno has p x timelength data
% in nominal conditions, the innovation vector should be zero mean white
% noise with covariance, V
% V = zeros( length(mrsInno(:,1)) ,1);
% for p = 1:length(mrsInno(:,1))
%     row_mean = mean(mrsInno(p,:));
%     V(p) = var(mrsInno(p,:));
% end
% 
% % My Chi-Squared Test (Siegert, 2016) 
% % for the whole time series
% df = p; % degrees of freedom
% S = diag(V'); % covariance
% Sinv = inv(S);
% test_stat = zeros(2,length(dateTvec));
% for i = 1:length(dateTvec)
%     y =  mrsInno(:,i);
%     test_stat(1,i) = y'*Sinv*y;
%     alpha = 0.001; % 95% confidence
% %     C(df,:) = chi2inv(1-[0.995 0.99 0.975 0.95 0.9 0.1 0.05 0.025 0.01 0.005],df);
%     X2(i) = chi2inv(1-alpha,df+4);
%     test_stat(2,i) = test_stat(1,i) > X2(i);
%     if mod(i,5000)==0
%         disp(i)
%     end
% end



% mrsInno has p x timelength data
% in nominal conditions, the innovation vector should be zero mean white
% noise with covariance, V
% V = zeros( length(mrsInno(:,1)) ,1);
% for m = 1:length(mrsInno(:,1))
%     row_mean = mean(mrsInno(m,:));
%     V(m) = var(mrsInno(m,:));
% end

% %% My Chi-Squared Test (Siegert, 2016) SIMPLE window
% % for the whole time series
% windowsize = 10000;
% df = m; % degrees of freedom
% test_stat = zeros(2,length(dateTvec));
% for i = windowsize+1:length(dateTvec)
%     S = diag(var(mrsInno(:,i-windowsize:i),0,2)');
%     y =  mrsInno(:,i);
%     test_stat(1,i) = y'*inv(S)*y;
%     alpha = 0.01; % 99% confidence
% %     C(df,:) = chi2inv(1-[0.995 0.99 0.975 0.95 0.9 0.1 0.05 0.025 0.01 0.005],df);
%     X2(i) = chi2inv(1-alpha,df);
%     test_stat(2,i) = test_stat(1,i) > X2(i);
%     if mod(i,5000)==0
%         disp(i)
%     end
% end

%% My Chi-Squared Test (Siegert, 2016) Eqn 33 window
% for the whole time series
windowsize = 10000;
df = length(mrsInno(:,1)); % degrees of freedom
test_stat = zeros(2,length(dateTvec));
for i = windowsize+1:length(dateTvec)
    S = diag(var(mrsInno(:,i-windowsize:i),0,2)');
    y =  mrsInno(:,i);
    test_stat(1,i) = y'*inv(S)*y;
    alpha = 0.01; % 99% confidence
    X2(i) = chi2inv(1-alpha,df);
    test_stat(2,i) = test_stat(1,i) > X2(i);
    if mod(i,5000)==0
        disp(i)
    end
end


%% plotting part
figure
yyaxis right
plot(dateTvec,test_stat(2,:),'*')
ylabel('IsRejected')

yyaxis left
plot(dateTvec,test_stat(1,:))
hold on
plot(dateTvec,X2)
ylabel('Test Statistic')

legend('Test Stat','Threshold','IsRejected','Location','best')


%% My GLRT (Siegert, 2016) 
% for the whole time series
windowsize = 10000;
df = length(mrsInno(:,1)); % degrees of freedom
test_stat = zeros(2,length(dateTvec));
for i = windowsize+1:length(dateTvec)
%     S = diag(var(mrsInno(:,i-windowsize:i),0,2)');
    S = C*sig_ekf(:,:,i)*C';
    y =  mrsInno(:,i);
    test_stat(1,i) = y'*inv(S)*y;
    alpha = 0.05; % 95% confidence
    X2(i) = chi2inv(1-alpha,df);
    test_stat(2,i) = test_stat(1,i) > X2(i);
    if mod(i,5000)==0
        disp(i)
    end
end







%% plotting part
figure
yyaxis right
plot(dateTvec,test_stat(2,:),'*')
ylabel('IsRejected')

yyaxis left
plot(dateTvec,test_stat(1,:))
hold on
plot(dateTvec,X2)
ylabel('Test Statistic')

legend('Test Stat','Threshold','IsRejected','Location','best')


