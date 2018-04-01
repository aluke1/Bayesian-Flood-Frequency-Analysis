clc
clear
%This program fits the stationary and non-stationary log-Pearson Type III
%distributions to peak streamflow records using the Bayesian methods
%presented in Luke, Adam, et al. "Predicting nonstationary flood frequencies: 
%Evidence supports an updated stationarity thesis in the United States." 
%Water Resources Research 53.7 (2017): 5469-5494.
%Code written by Adam Luke, November 2016: adluke61@gmail.com

%column 1 = water year, column 2 = annual max Q) 
data = load('record.txt');      %please remove zeros from observed Q
                                %no repeating water years in column 1
 
cf  = 0.3048^3;                 %conversion factor to convert input Q to m^3/s 
%cf = 1;
Mj  = 2;                        %1 for ST LPIII, 2 for NS LPIII with linear trend in mu
y_r = 0;                        %Regional estimate of gamma (skew coefficient)
SD_yr = 0.55;                   %Standard deviation of the regional estimate

%Prior distributions (input MATLAB abbreviation of distribution name used in 
%function call, i.e 'norm' for normal distribution as in 'normpdf')
marg_prior{1,1} = 'norm'; 
marg_prior{1,2} = 'unif'; 
marg_prior{1,3} = 'unif'; 
marg_prior{1,4} = 'unif';

%Hyper-parameters of prior distributions (input in order of use with 
%function call, i.e [mu, sigma] for normpdf(mu,sigma))
marg_par(:,1) = [y_r, SD_yr]';  %mean and std of informative prior on gamma 
marg_par(:,2) = [0, 2]';        %lower and upper bound of uniform prior on scale
marg_par(:,3) = [-10, 10]';     %lower and upper bound of uniform prior on location
marg_par(:,4) = [-0.15, 0.15]'; %lower and upper bound of uniform prior on trend 

%DREAM_(ZS) Variables
if Mj == 1; d = 3; end          %define number of parameters based on model
if Mj == 2; d = 4; end 
N = 3;                          %number of Markov chains 
T = 8000;                       %number of generations

%create function to initialize from prior
prior_draw = @(r,d)prior_rnd(r,d,marg_prior,marg_par); 
%create function to compute prior density 
prior_density = @(params)prior_pdf(params,d,marg_prior,marg_par);
%create function to compute unnormalized posterior density 
post_density = @(params)post_pdf(params,data,cf,Mj,prior_density);

%call the DREAM_ZS algorithm 
%Markov chains | post. density | archive of past states
[x,              p_x,          Z] = dream_zs(prior_draw,post_density,N,T,d,marg_prior,marg_par); 
%% Post Processing and Figures

%options:
%Which mu_t for calculating return level vs. return period? 
t = data(:,1) - min(data(:,1));                              %time (in years from the start of the fitting period)
idx_mu_n = size(t,1);                                        %calculates and plots RL vs RP for mu_t associated with t(idx_mu_n) 
                                                             %(idx_mu_n = size(t,1) for uST distribution) 
%Which return level for denisty plot? 
sRP = 100;                                                   %plots density of return level estimates for this return period

%Which return periods for output table?                      %outputs table of return level estimates for these return periods
RP_out =[200; 100; 50; 25; 10; 5; 2]; 
%end options 

%apply burn in (use only half of each chain) and rearrange chains to one sample 
x1 = x(round(T/2)+1:end,:,:);                                %burn in    
p_x1 = p_x(round(T/2)+1:end,:,:); 
post_sample = reshape(permute(x1,[1 3 2]),size(x1,1)*N,d,1); %columns are marginal posterior samples                                                           
sample_density = reshape(p_x1,size(p_x1,1)*N,1);             %corresponding unnormalized density 

%find MAP estimate of theta 
idx_theta_MAP = max(find(sample_density == max(sample_density))); 
theta_MAP = post_sample(idx_theta_MAP,:);                    %MAP parameter estimate 

%Compute mu as a function of time and credible intervals  
if Mj == 1; mu_t = repmat(post_sample(:,3),1,length(t));end  %ST model, mu is constant
if Mj == 2; mu_t = repmat(post_sample(:,3),1,length(t)) + post_sample(:,4)*t';end %NS mu = f(t)
MAP_mu_t = mu_t(idx_theta_MAP,:);                            %MAP estimate of the location parameter
low_mu_t = prctile(mu_t,2.5,1);                              %2.5 percentile of location parameter
high_mu_t = prctile(mu_t,97.5,1);                            %97.5 percentile of location parameter

%compute quantiles of the LPIII distribution 
p = 0.01:0.005:0.995;                                        %1st - 99.5th quantile (1 - 200 year RP)
a=1;
RLs = nan(size(post_sample,1),size(p,2));
for i = 1:size(post_sample,1);                               %compute return levels for each posterior sample
    RLs(i,:) = lp3inv(p,post_sample(i,1),post_sample(i,2),mu_t(i,idx_mu_n)); 
    a = a+1;
    if a == round(size(post_sample,1)/10) || i == size(post_sample,1);
    clc
    disp(['Calculating Return Levels ' num2str(round(i/size(post_sample,1)*100)) '% complete'])
    a = 1; 
    end
end
MAP_RL = RLs(idx_theta_MAP,:);                               %Return levels associated with MAP parameter estimate
low_RL = prctile(RLs,2.5,1);                                 %2.5 percentile of return level estimates
high_RL = prctile(RLs,97.5,1);                               %97.5 percentile of return level estimates
    
% Figure 1:  MAP and confidence interval of the location parameter vs data
figure
plot(data(:,1),data(:,2).*cf,'k-')
hold on 
plot(data(:,1),10.^MAP_mu_t,'r-') 
plot(data(:,1),10.^low_mu_t,'r--')
plot(data(idx_mu_n,1),10.^MAP_mu_t(1,idx_mu_n),'kx','linewidth',2)
plot(data(:,1),10.^high_mu_t,'r--')
l1 = legend('Annual Max. Q','MAP \mu_{\it{t}}','95% Credible Interval','MAP \mu_{\it{n}}'); 
l1.Location = 'Northwest';
xlim([min(data(:,1)) max(data(:,1))])
xlabel('Time (years)'); ylabel('Q (m^3/s)')

%Figue 2: Return level vs Return period plot for mu associated with t(idx_mu_n)
figure
plot(1./(1-p),MAP_RL,'k') 
hold on 
plot(1./(1-p),low_RL,'k--')  
plot(1./(1- [1:length(data)]./(length(data)+1)), sort(data(:,2)).*cf,'ro')
plot(1./(1-p),high_RL,'k--') 
xlabel('Return Period (years)'); ylabel('Return Level (m^3/s)'); 
l2 = legend('MAP Estimate','95% Credible Interval','Empirical (Stationary)');
l2.Location = 'northwest'; 
title(['Return Level vs Return Period (' num2str(data(idx_mu_n,1)) ')'])

%Figure 3: Distribution of selected return period 
idx_sRP = find(p == 1 - 1/sRP);
x_RL = min(RLs(:,idx_sRP)):max(RLs(:,idx_sRP));
y_RL = ksdensity(RLs(:,idx_sRP),x_RL);
figure
plot(x_RL,y_RL,'k');
xlabel([num2str(sRP) ' Year Return Level (' num2str(data(idx_mu_n,1)) ')']); ylabel('Density');

%output results to terminal 
clc
disp('MAP Parameter Estimate:')
if Mj == 1; disp(['skew = ' num2str(theta_MAP(1,1)) ' scale = '...
            num2str(theta_MAP(1,2)) ' location = ' num2str(theta_MAP(1,3))]);end
if Mj == 2; disp(['skew = ' num2str(theta_MAP(1,1)) ', scale = '...
            num2str(theta_MAP(1,2)) ', location = '     ...
            num2str(theta_MAP(1,3)) ', trend = ' num2str(theta_MAP(1,4))]);end
for i = 1:length(RP_out); 
idx_p= find(p == 1 - 1./RP_out(i)); 
table(i,1) = RP_out(i);       table(i,2) = MAP_RL(1,idx_p); 
table(i,3) = low_RL(1,idx_p); table(i,4) = high_RL(1,idx_p); 
end
disp(' ')
disp('           Return Level Estimates (cms)')     
disp('T(yrs)|  MAP est.    |  2.5 pctl.   |  97.5 pctl. ')
disp(num2str(table))









