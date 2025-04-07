tic

% @ Gerardo Chowell
% Last modification: 08/08/2017

%% 1) Fitting GGM model to early epidemic phase

clear

% Name of the data file containing the time series [ Time points| case incidence]
datafilename1='zika-daily-onset-antioquia.txt';

% load data file
data=load(datafilename1);

t_window=1:50; % the length of the early growth phase
disp(t_window)

alpha1=NaN; %0.4; % parameter for simple exponential smoothing; alpha1=NaN for unweighted NLSQ.

if isnan(alpha1) % weights are not assigned (unweighted NLSQ)
    weights=ones(length(t_window),1);
else
    weights=sesweights(length(t_window),alpha1); % weights according to simple exponential smoothing
end

% Parameters
r = 0.2;
p = 0.8;
a = 1;
K = 1000;
C0 = 20;

% Time vector
% (80 days)
timevect = 0:80;

% Generalized Richards Model (GRM) function
GRM = @(t, r, a, K, C0) K ./ (1 + ((K/C0 - 1) .* exp(-r * t)).^(1/a));

% Simulate data using GRM
simulated_data = GRM(timevect, r, a, K, C0);

% Add some noise to simulate real-world data (optional)
noise_level = 20;
simulated_data_noisy = simulated_data + noise_level * randn(size(simulated_data));
simulated_data_noisy = max(simulated_data_noisy, 0); % Ensure no negative cases

grmData = [timevect(:), simulated_data_noisy(:)];

data1=grmData(t_window, :); % select the early growth phase of the epidemic

DT=1; % temporal resolution of the data in days

flag1=2; % flag1=1 uses the exponential growth model; flag1=2 uses the generalized growth model

% create time vector in days
timevect=(data1(:,1))*DT;

I0=data1(1,2); % initial condition

z(1)=r;
z(2)=p;

LB=[0  0]; % lower bound for r and p
UB=[20  1]; % upper bound for r and p

% Optimization settings
options = optimoptions('lsqcurvefit', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'iter', ...
    'FiniteDifferenceType', 'central');

% Ensure sizes of inputs match
disp(size(timevect));
disp(size(weights.*data1(:,2)));

% Call lsqcurvefit
[P, resnorm, residual] = lsqcurvefit(@plotSimpleGrowth1, z, timevect, weights.*data1(:,2), LB, UB, options, I0, flag1, weights);

% P is the vector with the estimated parameters
r_hat=P(1)
p_hat=P(2)

% initial number of cases
IC(1)=I0

% obtain the best model fit
[~,F]=ode45(@simpleGrowth,timevect,IC,[],r_hat,p_hat,flag1);

% F contains the cumulative curve

% obtaian the incidence curve from the derivative of F
incidence1=[F(1,1);diff(F(:,1))];

% f contains incidence curve
f=incidence1;

% plot the model fit and data
figure
subplot(1,2,1)
line2=plot(timevect,incidence1,'r')
set(line2,'LineWidth',2)
hold on

line2=plot(grmData(:,1)*DT,grmData(:,2),'bo')
set(line2,'LineWidth',2)

xlabel('Time (days)');
ylabel('Incidence')
legend('Model','Data')

set(gca,'FontSize', 24);
set(gcf,'color','white')

% plot the residuals
subplot(1,2,2)
resid=(incidence1-data1(:,2));
scaledresid=resid./std(resid);

stem(timevect,resid)

xlabel('Time (days)');
ylabel('Residuals')

set(gca,'FontSize', 24);
set(gcf,'color','white')

% goodness of fit metrics
SSE=sum((incidence1-data1(:,2)).^2);


RMSE=sqrt(mean((incidence1-data1(:,2)).^2))
MAE=mean(abs(incidence1-data1(:,2)))
nonzero_idx = data1(:,2) ~= 0; % Indices where observed values are nonzero
MAPE = mean(abs((incidence1(nonzero_idx) - data1(nonzero_idx,2)) ./ data1(nonzero_idx,2))) * 100;

%% 2) Derive parameter uncertainty

% no weights used for generating parameter uncertainty
weights=ones(length(t_window),1);

% number of bootstrap realizations
M=100;

%distribution for data variation (dist1=1 for Poisson, dist1=2 for negative
%binomial, variance is mean1*factor1
dist1=1;
factor1=1;

% parameter estimates
Ptrue=[r_hat p_hat];

% F contains the cumulative number of cases derived from fitting the model to the data
% f contains incidence curve

Phatss=zeros(M,2); % vector that will store parameter estimates from each realization
SSEs=zeros(M,1); % SSE for each realization
curves=[]; 

% generate each of the M realizations
for realization=1:M
    
    % bootstrap realization is generated
    
    yirData=zeros(length(F),1);
    
    yirData(1)=F(1);
    
    
    if dist1==1         %Poisson dist. error structure
        for t=2:length(F)
            newcases_t=F(t)-F(t-1);
            yirData(t,1)=poissrnd(newcases_t,1,1);
        end
        
    elseif dist1==2         % Negative binomial dist. error structure
        for t=2:length(F)
            newcases_t=F(t)-F(t-1);
            mean1=newcases_t;
            var1=mean1*factor1
            
            param2=mean1/var1;
            param1=mean1*param2/(1-param2);
            
            yirData(t,1)=nbinrnd(param1,param2,1,1);
        end
        
    end
    
    % store realiation of the epidemic curve
    curves=[curves yirData];

    % Estimate parameters from bootstrap realization
    
    % initial values for parameters
    r=0.5;
    p=0.5;
   
    z(1)=r;% initial guesses
    z(2)=p;
    
    LB=[0  0]; % lower bound for beta and gamma
    UB=[20  1]; % upper bound for beta and gamma
    
    I0=data1(1,2); % initial condition
    
    hold on
        
    options = optimoptions('lsqcurvefit', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'iter', ...
    'FiniteDifferenceType', 'central');

% Ensure sizes of inputs match
disp(size(timevect));
disp(size(weights.*data1(:,2)));

% Call lsqcurvefit
[P, resnorm, residual] = lsqcurvefit(@plotSimpleGrowth1, z, timevect, weights.*yirData, LB, UB, options, I0, flag1, weights);
        
    % P is the vector with the estimated parameters from the realization
    r_hat=P(1);
    p_hat=P(2);
    
    % Stores estimated parameters for each realization of a total of M realizations
    Phatss(realization,:)=P; 
    
end

% estimate mean and 95% CI from distribution of parameter estimates
param_r=[mean(Phatss(:,1)) plims(Phatss(:,1),0.025) plims(Phatss(:,1),0.975)];
param_p=[mean(Phatss(:,2)) plims(Phatss(:,2),0.025) plims(Phatss(:,2),0.975)];

cad1=strcat('r=',num2str(param_r(end,1),2),' (95% CI:',num2str(param_r(end,2),2),',',num2str(param_r(end,3),2),')')
cad2=strcat('p=',num2str(param_p(end,1),2),' (95% CI:',num2str(param_p(end,2),2),',',num2str(param_p(end,3),2),')')

% plot empirical distributions of r and p
figure(100)
subplot(2,2,1)
hist(Phatss(:,1))
xlabel('r','fontweight','bold','fontsize',20);
ylabel('Frequency','fontweight','bold','fontsize',20)

title(cad1)

set(gca, 'FontSize', 24, 'LineWidth', 2, 'FontWeight', 'bold');
set(gcf,'color','white')

subplot(2,2,2)
hist(Phatss(:,2))
xlabel('p','fontweight','bold','fontsize',20);
ylabel('Frequency','fontweight','bold','fontsize',20)

title(cad2)

set(gca, 'FontSize', 24, 'LineWidth', 2, 'FontWeight', 'bold');
set(gcf,'color','white')

% plot model fit and uncertainty around model fit

param_r_time2=[];
param_p_time2=[];

subplot(2,2,[3 4])
plot(timevect,curves,'c-')
hold on

%line1=plot(timevect,mean(curves,2),'r-')
line1=plot(timevect,plims(curves',0.5),'r-')
set(line1,'LineWidth',2)

line1=plot(timevect,plims(curves',0.025),'r--')
set(line1,'LineWidth',2)

line1=plot(timevect,plims(curves',0.975),'r--')
set(line1,'LineWidth',2)

line1=plot(timevect,data1(:,2),'bo')
set(line1,'LineWidth',2)

hold on

xlabel('Time (days)');
ylabel('Case incidence')

set(gca,'FontSize', 24);
set(gcf,'color','white')
    
            
%% 3) Generate short-term forecasts



% forecast horizon
forecastingperiod=30; % number of data points ahead

timevect2=(0:t_window(end)-1+forecastingperiod)*DT;

% vector to store forecast curves
curvesforecasts1=[];

% generate forecast curves from each bootstrap realization
for realization=1:M
    
    r_hat=Phatss(realization,1);
    p_hat=Phatss(realization,2);
    
    [~,x]=ode45(@simpleGrowth,timevect2,IC,[],r_hat,p_hat,flag1);
    
    incidence1=[x(1,1);diff(x(:,1))];
    
    curvesforecasts1=[curvesforecasts1 incidence1];
    
end

% plot forecast curves
plot(timevect2,curvesforecasts1,'c')
hold on

line1=plot(timevect2,plims(curvesforecasts1',0.5),'r')
set(line1,'LineWidth',2)

line1=plot(timevect2,plims(curvesforecasts1',0.025),'r--')
set(line1,'LineWidth',2)

line1=plot(timevect2,plims(curvesforecasts1',0.975),'r--')
set(line1,'LineWidth',2)

xlabel('\fontsize{24}Time (days)');
ylabel('\fontsize{24}Case incidence')

% plot time series data
line1=plot(grmData(:,1)*DT,grmData(:,2),'ko')
set(line1,'LineWidth',2)

% plot vertical line separating calibration versus forecast periods
line2=[timevect(end) 0;timevect(end) max(grmData(:,2))*2];

%axis([timevect(1) timevect2(end) 0 max(grmData(:,2))*2])

line1=plot(line2(:,1),line2(:,2),'k--')
set(line1,'LineWidth',2)

set(gca,'FontSize',24)

% save data with results
save(strcat('Forecast-GGM-M-',num2str(M),'-tf-',num2str(t_window(end)),'-horizon-',num2str(forecastingperiod),'-flag-',num2str(flag1),'-dist-',num2str(dist1),'-factor-',num2str(factor1),'-w-',num2str(alpha1),'-',datafilename1(1:end-4),'.mat'),'-mat')

disp('Goodness of fit metrics:');
disp(['SSE: ', num2str(SSE)]);
disp(['RMSE: ', num2str(RMSE)]);
disp(['MAE: ', num2str(MAE)]);
disp(['MAPE: ', num2str(MAPE)]);

toc