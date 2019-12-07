%% SUNDELL AND SAYLOR (2017): AN INVERSE MODEL FOR UNMIXING DETRITAL GEOCHRONOLOGY AGE DISTRIBUTIONS; GEOCHEMISTRY, GEOPHYSICS, GEOSYSTEMS. %%
%% MODEL VERSION 1.1.02 %%

clear all % Clear all existing variables
close all % Close all figures
clc % Clear MATLAB Command Window

%% REQUIRED USER SPECIFICATIONS %%

trials = 10000; % Number of randomly generated sets of weights used to scale source distributions for comparison to mixed distribution
threshold_D = 1; % Set KS test D threshold, only change to conserve memory (i.e., if running > 100,000 trials and you know there are fits better than 0.2, set at 0.2 to retain nothing greater)
threshold_V = 1; % Basically the same as D threshold, but typically set a little higher since V = max(CDF1-CDF2) + max(CDF2-CDF1)
threshold_R2 = 0; % Same basically as D and V, except better fits are closer to 1 rather than 0
PDP_min = 0; % Minimum x (age) for constructing and plotting probability density plots
PDP_max = 3000; % Maximum x (age) for constructing and plotting probability density plots
PDP_step = 1; % Interval for constructing and plotting probability density plots, smaller steps make smoother PDPs, but cost compute memory
plotnum = 100; % number of best fits to retain (and plot) to calculate weight mean and standard dev. contribution from each source sample

%% OPTIONAL USER SPECIFICATIONS %%

% Random n ages resulting in N source distributions for each randomly generated set of weights during the Monte Carlo simulation. WARNING! Computationally intensive
Licht_approach = 0; % 0 = Do not use Licht et al. (2016) approach. 1 = Use approach during Monte Carlo model. 2 = Use approach following optimization (see below) 
Licht_N = 50; % N = the number of times each source distribution is subsampled (see Licht et al. (2016) for details, N = 200 in that paper)
Licht_n = 100; % n = the total number of ages sum of subsampled source ages; n must be <= minimum source sample size (see Licht et al. (2016) for details, n = 800 in that paper)

% Optimization options
Optimize = 0; % 0 = Do not run optimization routine. 1 = Run optimization routine. Optimization constrints set based on best-fit mean and standard dev. results from Monte Carlo model, see Sundell and Saylor (2017) for details
Max_best_fits = 5; % If using optimization routine #1, this is the number of best fits used to set bounds on source proportions in subsequent iterations.

%% INPUT SAMPLE DATA, CONSTRUCT SOURCE AND MIXED SAMPLE DISTRIBUTIONS, SET EMPTY VARIABLES %%

global data
global N
global x
global x_all
global cdf_sink 
global cdf_source
global pdp_sink 
global pdp_source

% Open browser window to find and load sample data; data need to be organized into age-uncertainty pairs (ages, uncertainties, ages, uncertianties, ...) with headers (sample names must start with a letter); first 2 columns are the mixed sample, all others are sources
[filename, pathname] = uigetfile({'*'},'File Selector');
fullpathname = strcat(pathname, filename);
[numbers, text, data_tmp] = xlsread(fullpathname);
data = cell2mat(data_tmp(2:end,:));
[dataR,dataC]=size(data);
N = (dataC/2)-1; % Number of source samples
tic
% Make probability density plots
x = PDP_min:PDP_step:PDP_max;
for i = 1:N+1;
m = data(:,i*2-1);
m = m(isfinite(m(:,1)),:);
s = data(:,i*2);
s = s(isfinite(s(:,1)),:);
f = zeros(length(m),length(x));
for j = 1:length(m);
f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*PDP_step);
end
pdps(:,i) = ((sum(f, 1))/length(m)).';
end
pdp_sink = pdps(:,1); % Mixed sample PDP
pdp_source = pdps(:,2:end); % All source sample PDPs

% Make cumulative distribution functions
for i = 1:N+1
ages(:,i) = (data(:,i*2-1));
end
x_all = sort(nonzeros(reshape(ages,[numel(ages),1])));
x_all(isnan(x_all(:,1)),:) = [];
x_all(x_all>4500) = []; % Remove ages older than Earth
x_all(x_all<0) = []; % Remove future ages
binEdges    =  [-inf ; x_all ; inf];
for i=1:N+1
x1 = ages(:,i);
x1(isnan(x1(:,1)),:) = [];
binCounts(:,i)  =  histc(nonzeros(x1), binEdges, 1);
clear x1
end
for i=1:N+1
sumCounts(:,i)  =  cumsum(binCounts(:,i))./sum(binCounts(:,i));
end
clear binCounts 
for i=1:N+1
CDF(:,i)  =  sumCounts(1:end-1, i);
end
clear sumCounts
cdf_sink = CDF(2:end,1); % Mixed sample CDF
cdf_source = CDF(2:end,2:N+1); % All source sample CDFs
clear CDF

% Set zero variables
fit_D = 0; % Keep model results if <= threshold_D
fit_V = 0; % Keep model results if <= threshold_V
fit_R2 = 0; % Keep model results if >= threshold_R2
count = 0; % Keep track of trial number during Monte Carlo model

% Set empty variables; if no model fits exist from Monte Carlo model due to threshold settings then these will remain empty ([])
Passed_D = [];
Passed_V = [];
Passed_R2 = [];
Passed_cdf_D = [];
Passed_cdf_V = [];
Passed_pdp_R2 = [];
min_D = [];
min_V = [];
max_R2 = [];
Results_D = [];
Results_V = [];
Results_R2 = [];

figure;
hold on
colours = colormap(jet((N)));
for i = 1:N;
plot(x_all,cdf_source(:,i),'color',colours((i),:),'linewidth',1.5);
grid on
end
q1=plot(x_all,cdf_sink, 'k', 'linewidth', 2);
legend([q1],{'Target'})
title('Source and Target CDFs')
xlabel('Age (Ma)')
ylabel('Cumulative probability');
ax.YAxis.Exponent = -2;
hold off

figure;
hold on
colours = colormap(jet((N)));
for i = 1:N;
plot(x,pdp_source(:,i),'color',colours((i),:),'linewidth',1.5);
grid on
end
q1=plot(x,pdp_sink, 'k', 'linewidth', 2);
legend([q1],{'Target'})
title('Source and Target PDPs')
xlabel('Age (Ma)')
ylabel('Relative probability');
ax.YAxis.Exponent = -2;
hold off

%% MONTE CARLO MODEL: RANDOMLY GENERATE SETS OF WEIGHTS TO SCALE SOURCE DISTRIBUTIONS FOR COMPARISON TO MIXED SAMPLE DISTRIBUTION %%

while count < trials % Run Monte Carlo until reach number of user specified trials

count = count + 1 % Current trial during Monte Carlo model

% Determine randomly generated set of weights for each trial
x3 = 1:1:N-1;
samples = x3(randperm(N-1));
x3 = 1:1:N;
m = 1/N;
sampled_y = 0;
sampled_y(1,N+1) = 1;
y_line = samples.*m;
for i = 1:N-1
tmp_min = max(sampled_y(1,1:samples(1,i)));    
tmp_max = min(nonzeros(sampled_y(1,samples(1,i)+1:end))); 
if tmp_min < y_line(1,i) 
sampled_y(1,samples(1,i)+1) = y_line(1,i) + rand*(tmp_max-y_line(1,i));
else
sampled_y(1,samples(1,i)+1) = tmp_min + rand*(tmp_max-tmp_min);
end
end
weights_tmp = diff(sampled_y);
weights = weights_tmp(randperm(length(weights_tmp))); % Randomly generated set of weights for each trial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF USER-SPECIFIED Licht_approach = 0 THEN RANDOMLY SAMPLED WEIGHTS WILL SCALE ENTIRE SOURCE DISTRIBUTIONS FOR COMPARISON TO MIXED SAMPLE DISTRIBUTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Licht_approach == 0

% Weight and sum all source distributions
cdf = sum(cdf_source.*(repmat(weights,length(x_all),1)),2);
pdp = sum(pdp_source.*(repmat(weights,length(x),1)),2);

% Compare randomly weighted sources with mixed sample
D = max([max(cdf_sink - cdf)],[max(cdf - cdf_sink)]);
V = max(cdf_sink - cdf) + max(cdf - cdf_sink);
R2 = ((sum((pdp_sink - mean(pdp_sink)).*(pdp - mean(pdp))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp - mean(pdp)).*(pdp - mean(pdp)))))))*...
	((sum((pdp_sink - mean(pdp_sink)).*(pdp - mean(pdp))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp - mean(pdp)).*(pdp - mean(pdp)))))));

D_std = 0;
V_std = 0;
R2_std = 0;

end % End if Licht_approach = 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF USER-SPECIFIED Licht_approach = 1 THEN Licht_N SOURCE DISTRIBUTIONS ARE CREATED BY RANDOMLY SAMPLING SOURCE AGES BASED ON trial WEIGHT FOR A TOTAL OF Licht_n AGES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Licht_approach == 1

% First determine how many ages correspond to each source based on rounding the randomly determined trial weight 
num_Ages = round(weights.*Licht_n); 

% Make sure all randomly sampled groups of ages sum to user specified Licht_n; if rounding results in total of ages more or fewer than Licht_n add the difference to the source with the fewest or subtract difference from source with the most
if sum(num_Ages) < Licht_n
dif = Licht_n - sum(num_Ages);
[num num_idx] = min(num_Ages);
num_Ages(1,num_idx) = num_Ages(1,num_idx) + dif;
end
if sum(num_Ages) > Licht_n
dif = sum(num_Ages) - Licht_n;
[num num_idx] = max(num_Ages);
num_Ages(1,num_idx) = num_Ages(1,num_idx) - dif;
end

% Randomly sample ages from each source age distribution Licht_N times based on num_Ages (above) and concatenate into pairs of columns (ages, uncertainties, ages, uncertianties, ...)
for k = 1:N
l_tmp = nonzeros(data(:,k*2+1));
l_tmp(isnan(l_tmp(:,1)),:) = [];
l_ages(k,1) = length(l_tmp);
end
	if min(l_ages) < Licht_n
	err_dlg=errordlg('The number of subsamples must be less than or equal to smallest source sample because the model samples without replacement.','Hold on a sec...');
	waitfor(err_dlg);
	else
	end
rand_Ages = zeros(Licht_n,N*2,Licht_N);
rand_ages_tmp1 = [];
for p = 1:Licht_N
for i = 1:N
if num_Ages(1,i) > 0
data_tmp = data(:,2*i+1:2*i+2);
data_tmp = data_tmp(any(~isnan(data_tmp),2),:); % remove NaNs
data_tmp(all(data_tmp==0,2),:)=[];
sz = size(rand_ages_tmp1);
rand_ages_tmp1(sz(1,1)+1:sz(1,1)+num_Ages(1,i),p*2-1:p*2) = datasample(data_tmp,num_Ages(1,i),1, 'Replace', false);
end
end
end
for p = 1:Licht_N
rand_ages_tmp2 = rand_ages_tmp1(:,p*2-1:p*2);
rand_ages_tmp2(all(rand_ages_tmp2==0,2),:)=[];
rand_ages_all(:,p*2-1:p*2) = rand_ages_tmp2;
clear rand_ages_tmp2
end

% Combine mixed sample distribution and randomly generated model distribution to generate CDF curves with equally binned x axes (required to calculate D and V)
ages_licht = zeros(max([length(data(:,1)),Licht_n]),Licht_N+1);
ages_licht(1:length(data(:,1)),1) = data(:,1);
for p = 1:Licht_N
ages_licht(1:Licht_n,p+1) = rand_ages_all(:,p*2-1);
end
x_all = sort(nonzeros(reshape(ages_licht,[numel(ages_licht),1])));
x_all(isnan(x_all(:,1)),:) = [];
x_all(x_all>4500) = []; % Remove ages older than Earth
x_all(x_all<0) = []; % Remove future ages
binEdges    =  [-inf ; x_all ; inf];
for i=1:Licht_N+1
x1 = ages_licht(:,i);
x1(isnan(x1(:,1)),:) = [];
binCounts(:,i)  =  histc(nonzeros(x1), binEdges, 1);
clear x1
end
for i=1:Licht_N+1 
sumCounts(:,i)  =  cumsum(binCounts(:,i))./sum(binCounts(:,i));
end
for i=1:Licht_N+1
CDF(:,i)  =  sumCounts(1:end-1, i);
end
cdf_sink = CDF(2:end,1);
cdf_source = CDF(2:end,2:end);
cdf = mean(cdf_source,2);

% Generate source PDPs based on randomly sampled ages
for p = 1:Licht_N
m = rand_ages_all(:,p*2-1);
s = rand_ages_all(:,p*2);
pdp_rand_ages_all = ((sum(bsxfun(@times,1./(s.*sqrt(2*pi)),exp(bsxfun(@rdivide, -(bsxfun(@minus, x, m).^2), 2*s.^2))), 1)/length(m)).').*PDP_step;
pdp_source(:,p) = pdp_rand_ages_all;
end
pdp = mean(pdp_source,2);




% Compare each randomly generated age distribution to the mixed sample distribution for each trial weight
for p = 1:Licht_N
D_tmp(:,p) = max([max(cdf_sink - cdf_source(:,p))],[max(cdf_source(:,p) - cdf_sink)]);
V_tmp(:,p) = max(cdf_sink - cdf_source(:,p)) + max(cdf_source(:,p) - cdf_sink);
R2_tmp(:,p) = ((sum((pdp_sink - mean(pdp_sink)).*(pdp_source(:,p) - mean(pdp_source(:,p)))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp_source(:,p) - ...
	mean(pdp_source(:,p))).*(pdp_source(:,p) - mean(pdp_source(:,p))))))))*((sum((pdp_sink - mean(pdp_sink)).*(pdp_source(:,p) - mean(pdp_source(:,p)))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*...
	(pdp_sink - mean(pdp_sink))))*(sum((pdp_source(:,p) - mean(pdp_source(:,p))).*(pdp_source(:,p) - mean(pdp_source(:,p))))))));
end

% Calculate mean and stdev for each group of comparisons (total of Licht_N) within each trial; will be filtered below based on mean value

D = min(D_tmp);
V = min(V_tmp);
R2 = max(R2_tmp);

D_SE(count,1) = std(D_tmp)./sqrt(Licht_N);
V_SE(count,1) = std(V_tmp)./sqrt(Licht_N);
R2_SE(count,1) = std(R2_tmp)./sqrt(Licht_N);

D_std = std(D_tmp);
V_std = std(V_tmp);
R2_std = std(R2_tmp);

end % End if Licht_approach = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTER AND SORT MODEL RESULTS BASED ON USER SPECIFIED THRESHOLDS, DETERMINE BEST FIT SOURCE CONTRIBUTIONS, PLOT RESULTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Keep results if they pass user specified thresholds
if D <= threshold_D
fit_D = fit_D + 1;
Passed_cdf_D(:,fit_D) = cdf;
Passed_weights_D(fit_D, :) = weights;
Passed_D(fit_D, 1) = D;
Passed_x_all_D(:,fit_D) = x_all;
Passed_D_std(fit_D, 1) = D_std;
end

if V <= threshold_V
fit_V = fit_V + 1;
Passed_cdf_V(:,fit_V) = cdf;
Passed_weights_V(fit_V, :) = weights;
Passed_V(fit_V, 1) = V;
Passed_x_all_V(:,fit_V) = x_all;
Passed_V_std(fit_V, 1) = V_std;
end

if R2 >= threshold_R2
fit_R2 = fit_R2 + 1;
Passed_pdp_R2(:,fit_R2) = pdp;
Passed_weights_R2(fit_R2, :) = weights;
Passed_R2(fit_R2, 1) = R2;
Passed_R2_std(fit_R2,:) = R2_std;
end

end % End once reached number of user specified trials

% Sort and trim results to retain 'plotnum' number of best model fits
if length(Passed_D) >= plotnum
Results_D = sortrows([Passed_D, Passed_D_std, Passed_weights_D, transpose(Passed_cdf_D), transpose(Passed_x_all_D)], 1);
Results_D = Results_D(1:plotnum,:);
Results_D_mean_std = transpose([mean(Results_D(1:plotnum,3:N+2));std(Results_D(1:plotnum,3:N+2))]);
end

if length(Passed_V) >= plotnum
Results_V = sortrows([Passed_V, Passed_V_std, Passed_weights_V, transpose(Passed_cdf_V), transpose(Passed_x_all_V)], 1);
Results_V = Results_V(1:plotnum,:);
Results_V_mean_std = transpose([mean(Results_V(1:plotnum,3:N+2));std(Results_V(1:plotnum,3:N+2))]);
end

if length(Passed_R2) >= plotnum
Results_R2 = flipud(sortrows([Passed_R2, Passed_R2_std, Passed_weights_R2, transpose(Passed_pdp_R2)], 1));
Results_R2 = Results_R2(1:plotnum,:);
Results_R2_mean_std = transpose([mean(Results_R2(1:plotnum,3:N+2));std(Results_R2(1:plotnum,3:N+2))]);
end

% Concatenate mean and standard dev. of best fit source contributions for all metrics
if length(Passed_D) >= plotnum && length(Passed_V) >= plotnum && length(Passed_R2) >= plotnum
Results_All_mean_std = [Results_D_mean_std, Results_V_mean_std, Results_R2_mean_std];
end

D_mean = mean(Results_D(:,1));
V_mean = mean(Results_V(:,1));
R2_mean = mean(Results_R2(:,1));

if Licht_approach == 0
D_std = std(Results_D(:,1));
V_std = std(Results_V(:,1));
R2_std = std(Results_R2(:,1));
else
for i = 1:length(Results_D(:,2))
D_err_prop(i,1) = Results_D(i,2)^2;
D_std = sqrt(sum(D_err_prop));
end
for i = 1:length(Results_V(:,2))
V_err_prop(i,1) = Results_V(i,2)^2;
V_std = sqrt(sum(V_err_prop));
end
for i = 1:length(Results_R2(:,2))
R2_err_prop(i,1) = Results_R2(i,2)^2;
R2_std = sqrt(sum(R2_err_prop));
end
end

% Calculate mean and stdev. of best fits and best overall fit. Print results in MATLAB command window
mean_stdev_MONTECARLO_D_V_R2 = [[D_mean; D_std], [V_mean; V_std], [R2_mean; R2_std]]
best_MONTECARLO_D_V_R2 = [min(Passed_D), min(Passed_V), max(Passed_R2)]
best_MONTECARLO_weights = [Results_D(1,3:N+2); Results_V(1,3:N+2); Results_R2(1,3:N+2)]

% Plot 'plotnum' number of best model fits for KS test D value
if length(Results_D(:,1)) >= plotnum
F1 = figure;
for i = 1:plotnum
plot(Results_D(i,N+3+length(Passed_cdf_D(:,1)):end),Results_D(i,N+3:N+2+length(Passed_cdf_D(:,1))), 'color', [0.01 .68 .3] , 'linewidth', 1)
hold on
end
cdf_sink_ages = nonzeros(data(:,1));
cdf_sink_ages(isnan(cdf_sink_ages(:,1)),:) = [];
j=cdfplot(cdf_sink_ages(:,1));
set(j,'color','k','linewidth',2)
hold off
axis([PDP_min PDP_max 0 1])
title('KS test D Statistic')
end

% Plot 'plotnum' number of best model fits for Kuiper test V value
if length(Results_V(:,1)) >= plotnum
F2 = figure;
for i = 1:plotnum
plot(Results_V(i,N+3+length(Passed_cdf_V(:,1)):end),Results_V(i,N+3:N+2+length(Passed_cdf_V(:,1))), 'color', [1 .1 0], 'linewidth', 1)
hold on
end
cdf_sink_ages = nonzeros(data(:,1));
cdf_sink_ages(isnan(cdf_sink_ages(:,1)),:) = [];
j=cdfplot(cdf_sink_ages(:,1));
set(j,'color','k','linewidth',2)
hold off
axis([PDP_min PDP_max 0 1])
title('Kuier test V Statistic')
end

% Plot 'plotnum' number of best model fits for cross correlation (R2) of PDPs
if length(Results_R2) >= plotnum
F3 = figure;
for i = 1:plotnum
plot(x,Results_R2(i,N+3:end), 'color', [.04 .31 1], 'linewidth', 1)
hold on
end
plot(x,pdp_sink, 'k', 'linewidth', 2)
hold off
axis([PDP_min PDP_max 0 max(pdp_sink)+max(pdp_sink)*.25])
title('R2 CDF weighting')
end


% Visualize Monte Carlo plots (remove curly brackets to plot)

%{

KS_D = figure;
colours = colormap(jet((N)));
colorbar;
hold on
for j = 1:length(Passed_D)
for i = 1:N
line(Passed_weights_D(j,i),Passed_D(j,1),'LineStyle','none','Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
	'MarkerSize',6);
end
end
figure(KS_D)
set(0, 'CurrentFigure', KS_D)
set(KS_D,'Visible','on')
ylabel('KS test D value');
xlabel('Source contribution');
title('Monte Carlo mixture modeling using KS test D value');
set(gca,'Ydir','reverse')
hold off

Kuiper_V = figure;
colours = colormap(jet((N)));
colorbar;
hold on
for j = 1:length(Passed_V)
for i = 1:N
line(Passed_weights_V(j,i),Passed_V(j,1),'LineStyle','none','Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
	'MarkerSize',6);
end
end
figure(Kuiper_V)
set(0, 'CurrentFigure', Kuiper_V)
set(Kuiper_V,'Visible','on')
ylabel('Kuiper test V value');
xlabel('Source contribution');
title('Monte Carlo mixture modeling using Kuiper test V value');
set(gca,'Ydir','reverse')
hold off

R2 = figure;
colours = colormap(jet((N)));
colorbar;
hold on
for j = 1:length(Passed_R2)
for i = 1:N
line(Passed_weights_R2(j,i),Passed_R2(j,1),'LineStyle','none','Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
	'MarkerSize',6);
end
end
figure(R2)
set(R2,'Visible','on')
ylabel('Cross-correlation coefficient');
xlabel('Source contribution');
title('Monte Carlo mixture modeling using Cross-correlation coefficient');
hold off

%}

%% OPTIMIZATION OPTION #1, CONSTRAINTS SET BY MEAN AND STANDARD DEVIATION OF BEST FIT MONTE CARLO MODEL RESULTS %%

if Optimize == 1

x = PDP_min:PDP_step:PDP_max;
for i = 1:N+1;
m = data(:,i*2-1);
m = m(isfinite(m(:,1)),:);
s = data(:,i*2);
s = s(isfinite(s(:,1)),:);
f = zeros(length(m),length(x));
for j = 1:length(m);
f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*PDP_step);
end
pdps(:,i) = ((sum(f, 1))/length(m)).';
end
pdp_sink = pdps(:,1); % Mixed sample PDP
pdp_source = pdps(:,2:end); % All source sample PDPs

% Make cumulative distribution functions
binCounts = [];
sumCounts = [];
sumCounts = [];
CDF = [];
for i = 1:N+1
ages(:,i) = (data(:,i*2-1));
end
x_all = sort(nonzeros(reshape(ages,[numel(ages),1])));
x_all(isnan(x_all(:,1)),:) = [];
x_all(x_all>4500) = []; % Remove ages older than Earth
x_all(x_all<0) = []; % Remove future ages
binEdges    =  [-inf ; x_all ; inf];
for i=1:N+1
x1 = ages(:,i);
x1(isnan(x1(:,1)),:) = [];
binCounts(:,i)  =  histc(nonzeros(x1), binEdges, 1);
clear x1
end
for i=1:N+1
sumCounts(:,i)  =  cumsum(binCounts(:,i))./sum(binCounts(:,i));
end
clear binCounts 
for i=1:N+1
CDF(:,i)  =  sumCounts(1:end-1, i);
end
clear sumCounts
cdf_sink = CDF(2:end,1); % Mixed sample CDF
cdf_source = CDF(2:end,2:N+1); % All source sample CDFs
clear CDF

stdev = 1; 
iteration_int = [10, 5, 2, 1]; % Grid spacing (in percent). First find best fits at 10% spacing, then 5%, then 2%, then 1%.

for t = 1:length(iteration_int)

INT = iteration_int(1,t); % Set grid spacing

nodes = 0; % Reset number of different weight possibilities

if t == 1

% Set initial grid spacing based on Monte Carlo mean and stdev. and force to be divisibl by 10. For example, if mean weight is 22+/-5 (range of 17 to 27), grid for that source would be set at 10, 20, 30 for iteration 1.
a_min_D = round((Results_D_mean_std(1,1) - Results_D_mean_std(1,2)*stdev)*100);
a_max_D = round((Results_D_mean_std(1,1) + Results_D_mean_std(1,2)*stdev)*100);
a_min_D_rem = rem(a_min_D(1,1),10);
a_max_D_rem = 10 - rem(a_max_D(1,1),10);
if a_min_D < 0
a_min_D = 0;
a_min_D_rem = 0;
end
if a_max_D >= 100
a_max_D = 100;
a_max_D_rem = 0;
end
a_D = (a_min_D-a_min_D_rem):INT:(a_max_D+a_max_D_rem);

b_min_D = round((Results_D_mean_std(2,1) - Results_D_mean_std(2,2)*stdev)*100);
b_max_D = round((Results_D_mean_std(2,1) + Results_D_mean_std(2,2)*stdev)*100);
b_min_D_rem = rem(b_min_D(1,1),10);
b_max_D_rem = 10 - rem(b_max_D(1,1),10);
if b_min_D < 0
b_min_D = 0;
b_min_D_rem = 0;
end
if b_max_D >= 100
b_max_D = 100;
b_max_D_rem = 0;
end
b_D = (b_min_D-b_min_D_rem):INT:(b_max_D+b_max_D_rem);

if N == 2
weights_allcomb_D = allcomb(a_D,b_D);
end

if N >= 3 
c_min_D = round((Results_D_mean_std(3,1) - Results_D_mean_std(3,2)*stdev)*100);
c_max_D = round((Results_D_mean_std(3,1) + Results_D_mean_std(3,2)*stdev)*100);
c_min_D_rem = rem(c_min_D(1,1),10);
c_max_D_rem = 10 - rem(c_max_D(1,1),10);
if c_min_D < 0
c_min_D = 0;
c_min_D_rem = 0;
end
if c_max_D >= 100
c_max_D = 100;
c_max_D_rem = 0;
end
c_D = (c_min_D-c_min_D_rem):INT:(c_max_D+c_max_D_rem);
end
if N == 3
weights_allcomb_D = allcomb(a_D,b_D,c_D);
end

if N >= 4
d_min_D = round((Results_D_mean_std(4,1) - Results_D_mean_std(4,2)*stdev)*100);
d_max_D = round((Results_D_mean_std(4,1) + Results_D_mean_std(4,2)*stdev)*100);
d_min_D_rem = rem(d_min_D(1,1),10);
d_max_D_rem = 10 - rem(d_max_D(1,1),10);
if d_min_D < 0
d_min_D = 0;
d_min_D_rem = 0;
end
if d_max_D >= 100
d_max_D = 100;
d_max_D_rem = 0;
end
d_D = (d_min_D-d_min_D_rem):INT:(d_max_D+d_max_D_rem);
end
if N == 4
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D);
end

if N >= 5
e_min_D = round((Results_D_mean_std(5,1) - Results_D_mean_std(5,2)*stdev)*100);
e_max_D = round((Results_D_mean_std(5,1) + Results_D_mean_std(5,2)*stdev)*100);
e_min_D_rem = rem(e_min_D(1,1),10);
e_max_D_rem = 10 - rem(e_max_D(1,1),10);
if e_min_D < 0
e_min_D = 0;
e_min_D_rem = 0;
end
if e_max_D >= 100
e_max_D = 100;
e_max_D_rem = 0;
end
e_D = (e_min_D-e_min_D_rem):INT:(e_max_D+e_max_D_rem);
end
if N == 5
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D);
end

if N >= 6
f_min_D = round((Results_D_mean_std(6,1) - Results_D_mean_std(6,2)*stdev)*100);
f_max_D = round((Results_D_mean_std(6,1) + Results_D_mean_std(6,2)*stdev)*100);
f_min_D_rem = rem(f_min_D(1,1),10);
f_max_D_rem = 10 - rem(f_max_D(1,1),10);
if f_min_D < 0
f_min_D = 0;
f_min_D_rem = 0;
end
if f_max_D >= 100
f_max_D = 100;
f_max_D_rem = 0;
end
f_D = (f_min_D-f_min_D_rem):INT:(f_max_D+f_max_D_rem);
end
if N == 6
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D);
end

if N >= 7
g_min_D = round((Results_D_mean_std(7,1) - Results_D_mean_std(7,2)*stdev)*100);
g_max_D = round((Results_D_mean_std(7,1) + Results_D_mean_std(7,2)*stdev)*100);
g_min_D_rem = rem(g_min_D(1,1),10);
g_max_D_rem = 10 - rem(g_max_D(1,1),10);
if g_min_D < 0
g_min_D = 0;
g_min_D_rem = 0;
end
if g_max_D >= 100
g_max_D = 100;
g_max_D_rem = 0;
end
g_D = (g_min_D-g_min_D_rem):INT:(g_max_D+g_max_D_rem);
end
if N == 7
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D);
end

if N >= 8
h_min_D = round((Results_D_mean_std(8,1) - Results_D_mean_std(8,2)*stdev)*100);
h_max_D = round((Results_D_mean_std(8,1) + Results_D_mean_std(8,2)*stdev)*100);
h_min_D_rem = rem(h_min_D(1,1),10);
h_max_D_rem = 10 - rem(h_max_D(1,1),10);
if h_min_D < 0
h_min_D = 0;
h_min_D_rem = 0;
end
if h_max_D >= 100
h_max_D = 100;
h_max_D_rem = 0;
end
h_D = (h_min_D-h_min_D_rem):INT:(h_max_D+h_max_D_rem);
end
if N == 8
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D);
end

if N >= 9
l_min_D = round((Results_D_mean_std(9,1) - Results_D_mean_std(9,2)*stdev)*100);
l_max_D = round((Results_D_mean_std(9,1) + Results_D_mean_std(9,2)*stdev)*100);
l_min_D_rem = rem(l_min_D(1,1),10);
l_max_D_rem = 10 - rem(l_max_D(1,1),10);
if l_min_D < 0
l_min_D = 0;
l_min_D_rem = 0;
end
if l_max_D >= 100
l_max_D = 100;
l_max_D_rem = 0;
end
l_D = (l_min_D-l_min_D_rem):INT:(l_max_D+l_max_D_rem);
end
if N == 9
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D,l_D);
end

if N >= 10
q_min_D = round((Results_D_mean_std(10,1) - Results_D_mean_std(10,2)*stdev)*100);
q_max_D = round((Results_D_mean_std(10,1) + Results_D_mean_std(10,2)*stdev)*100);
q_min_D_rem = rem(q_min_D(1,1),10);
q_max_D_rem = 10 - rem(q_max_D(1,1),10);
if q_min_D < 0
q_min_D = 0;
q_min_D_rem = 0;
end
if q_max_D >= 100
q_max_D = 100;
q_max_D_rem = 0;
end
q_D = (q_min_D-q_min_D_rem):INT:(q_max_D+q_max_D_rem);
end
if N == 10
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D,l_D,q_D);
end

weights_D = weights_allcomb_D;
number_of_weights_D=length(weights_allcomb_D);

a_min_V = round((Results_V_mean_std(1,1) - Results_V_mean_std(1,2)*stdev)*100);
a_max_V = round((Results_V_mean_std(1,1) + Results_V_mean_std(1,2)*stdev)*100);
a_min_V_rem = rem(a_min_V(1,1),10);
a_max_V_rem = 10 - rem(a_max_V(1,1),10);
if a_min_V < 0
a_min_V = 0;
a_min_V_rem = 0;
end
if a_max_V >= 100
a_max_V = 100;
a_max_V_rem = 0;
end
a_V = (a_min_V-a_min_V_rem):INT:(a_max_V+a_max_V_rem);

b_min_V = round((Results_V_mean_std(2,1) - Results_V_mean_std(2,2)*stdev)*100);
b_max_V = round((Results_V_mean_std(2,1) + Results_V_mean_std(2,2)*stdev)*100);
b_min_V_rem = rem(b_min_V(1,1),10);
b_max_V_rem = 10 - rem(b_max_V(1,1),10);
if b_min_V < 0
b_min_V = 0;
b_min_V_rem = 0;
end
if b_max_V >= 100
b_max_V = 100;
b_max_V_rem = 0;
end
b_V = (b_min_V-b_min_V_rem):INT:(b_max_V+b_max_V_rem);

if N == 2
weights_allcomb_V = allcomb(a_V,b_V);
end

if N >= 3 
c_min_V = round((Results_V_mean_std(3,1) - Results_V_mean_std(3,2)*stdev)*100);
c_max_V = round((Results_V_mean_std(3,1) + Results_V_mean_std(3,2)*stdev)*100);
c_min_V_rem = rem(c_min_V(1,1),10);
c_max_V_rem = 10 - rem(c_max_V(1,1),10);
if c_min_V < 0
c_min_V = 0;
c_min_V_rem = 0;
end
if c_max_V >= 100
c_max_V = 100;
c_max_V_rem = 0;
end
c_V = (c_min_V-c_min_V_rem):INT:(c_max_V+c_max_V_rem);
end
if N == 3
weights_allcomb_V = allcomb(a_V,b_V,c_V);
end

if N >= 4
d_min_V = round((Results_V_mean_std(4,1) - Results_V_mean_std(4,2)*stdev)*100);
d_max_V = round((Results_V_mean_std(4,1) + Results_V_mean_std(4,2)*stdev)*100);
d_min_V_rem = rem(d_min_V(1,1),10);
d_max_V_rem = 10 - rem(d_max_V(1,1),10);
if d_min_V < 0
d_min_V = 0;
d_min_V_rem = 0;
end
if d_max_V >= 100
d_max_V = 100;
d_max_V_rem = 0;
end
d_V = (d_min_V-d_min_V_rem):INT:(d_max_V+d_max_V_rem);
end
if N == 4
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V);
end

if N >= 5
e_min_V = round((Results_V_mean_std(5,1) - Results_V_mean_std(5,2)*stdev)*100);
e_max_V = round((Results_V_mean_std(5,1) + Results_V_mean_std(5,2)*stdev)*100);
e_min_V_rem = rem(e_min_V(1,1),10);
e_max_V_rem = 10 - rem(e_max_V(1,1),10);
if e_min_V < 0
e_min_V = 0;
e_min_V_rem = 0;
end
if e_max_V >= 100
e_max_V = 100;
e_max_V_rem = 0;
end
e_V = (e_min_V-e_min_V_rem):INT:(e_max_V+e_max_V_rem);
end
if N == 5
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V);
end

if N >= 6
f_min_V = round((Results_V_mean_std(6,1) - Results_V_mean_std(6,2)*stdev)*100);
f_max_V = round((Results_V_mean_std(6,1) + Results_V_mean_std(6,2)*stdev)*100);
f_min_V_rem = rem(f_min_V(1,1),10);
f_max_V_rem = 10 - rem(f_max_V(1,1),10);
if f_min_V < 0
f_min_V = 0;
f_min_V_rem = 0;
end
if f_max_V >= 100
f_max_V = 100;
f_max_V_rem = 0;
end
f_V = (f_min_V-f_min_V_rem):INT:(f_max_V+f_max_V_rem);
end
if N == 6
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V);
end

if N >= 7
g_min_V = round((Results_V_mean_std(7,1) - Results_V_mean_std(7,2)*stdev)*100);
g_max_V = round((Results_V_mean_std(7,1) + Results_V_mean_std(7,2)*stdev)*100);
g_min_V_rem = rem(g_min_V(1,1),10);
g_max_V_rem = 10 - rem(g_max_V(1,1),10);
if g_min_V < 0
g_min_V = 0;
g_min_V_rem = 0;
end
if g_max_V >= 100
g_max_V = 100;
g_max_V_rem = 0;
end
g_V = (g_min_V-g_min_V_rem):INT:(g_max_V+g_max_V_rem);
end
if N == 7
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V);
end

if N >= 8
h_min_V = round((Results_V_mean_std(8,1) - Results_V_mean_std(8,2)*stdev)*100);
h_max_V = round((Results_V_mean_std(8,1) + Results_V_mean_std(8,2)*stdev)*100);
h_min_V_rem = rem(h_min_V(1,1),10);
h_max_V_rem = 10 - rem(h_max_V(1,1),10);
if h_min_V < 0
h_min_V = 0;
h_min_V_rem = 0;
end
if h_max_V >= 100
h_max_V = 100;
h_max_V_rem = 0;
end
h_V = (h_min_V-h_min_V_rem):INT:(h_max_V+h_max_V_rem);
end
if N == 8
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V);
end

if N >= 9
l_min_V = round((Results_V_mean_std(9,1) - Results_V_mean_std(9,2)*stdev)*100);
l_max_V = round((Results_V_mean_std(9,1) + Results_V_mean_std(9,2)*stdev)*100);
l_min_V_rem = rem(l_min_V(1,1),10);
l_max_V_rem = 10 - rem(l_max_V(1,1),10);
if l_min_V < 0
l_min_V = 0;
l_min_V_rem = 0;
end
if l_max_V >= 100
l_max_V = 100;
l_max_V_rem = 0;
end
l_V = (l_min_V-l_min_V_rem):INT:(l_max_V+l_max_V_rem);
end
if N == 9
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V,l_V);
end

if N >= 10
q_min_V = round((Results_V_mean_std(10,1) - Results_V_mean_std(10,2)*stdev)*100);
q_max_V = round((Results_V_mean_std(10,1) + Results_V_mean_std(10,2)*stdev)*100);
q_min_V_rem = rem(q_min_V(1,1),10);
q_max_V_rem = 10 - rem(q_max_V(1,1),10);
if q_min_V < 0
q_min_V = 0;
q_min_V_rem = 0;
end
if q_max_V >= 100
q_max_V = 100;
q_max_V_rem = 0;
end
q_V = (q_min_V-q_min_V_rem):INT:(q_max_V+q_max_V_rem);
end
if N == 10
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V,l_V,q_V);
end

weights_V = weights_allcomb_V;
number_of_weights_V=length(weights_allcomb_V);

a_min_R2 = round((Results_R2_mean_std(1,1) - Results_R2_mean_std(1,2)*stdev)*100);
a_max_R2 = round((Results_R2_mean_std(1,1) + Results_R2_mean_std(1,2)*stdev)*100);
a_min_R2_rem = rem(a_min_R2(1,1),10);
a_max_R2_rem = 10 - rem(a_max_R2(1,1),10);
if a_min_R2 < 0
a_min_R2 = 0;
a_min_R2_rem = 0;
end
if a_max_R2 >= 100
a_max_R2 = 100;
a_max_R2_rem = 0;
end
a_R2 = (a_min_R2-a_min_R2_rem):INT:(a_max_R2+a_max_R2_rem);

b_min_R2 = round((Results_R2_mean_std(2,1) - Results_R2_mean_std(2,2)*stdev)*100);
b_max_R2 = round((Results_R2_mean_std(2,1) + Results_R2_mean_std(2,2)*stdev)*100);
b_min_R2_rem = rem(b_min_R2(1,1),10);
b_max_R2_rem = 10 - rem(b_max_R2(1,1),10);
if b_min_R2 < 0
b_min_R2 = 0;
b_min_R2_rem = 0;
end
if b_max_R2 >= 100
b_max_R2 = 100;
b_max_R2_rem = 0;
end
b_R2 = (b_min_R2-b_min_R2_rem):INT:(b_max_R2+b_max_R2_rem);

if N == 2
weights_allcomb_R2 = allcomb(a_R2,b_R2);
end

if N >= 3 
c_min_R2 = round((Results_R2_mean_std(3,1) - Results_R2_mean_std(3,2)*stdev)*100);
c_max_R2 = round((Results_R2_mean_std(3,1) + Results_R2_mean_std(3,2)*stdev)*100);
c_min_R2_rem = rem(c_min_R2(1,1),10);
c_max_R2_rem = 10 - rem(c_max_R2(1,1),10);
if c_min_R2 < 0
c_min_R2 = 0;
c_min_R2_rem = 0;
end
if c_max_R2 >= 100
c_max_R2 = 100;
c_max_R2_rem = 0;
end
c_R2 = (c_min_R2-c_min_R2_rem):INT:(c_max_R2+c_max_R2_rem);
end
if N == 3
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2);
end

if N >= 4
d_min_R2 = round((Results_R2_mean_std(4,1) - Results_R2_mean_std(4,2)*stdev)*100);
d_max_R2 = round((Results_R2_mean_std(4,1) + Results_R2_mean_std(4,2)*stdev)*100);
d_min_R2_rem = rem(d_min_R2(1,1),10);
d_max_R2_rem = 10 - rem(d_max_R2(1,1),10);
if d_min_R2 < 0
d_min_R2 = 0;
d_min_R2_rem = 0;
end
if d_max_R2 >= 100
d_max_R2 = 100;
d_max_R2_rem = 0;
end
d_R2 = (d_min_R2-d_min_R2_rem):INT:(d_max_R2+d_max_R2_rem);
end
if N == 4
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2);
end

if N >= 5
e_min_R2 = round((Results_R2_mean_std(5,1) - Results_R2_mean_std(5,2)*stdev)*100);
e_max_R2 = round((Results_R2_mean_std(5,1) + Results_R2_mean_std(5,2)*stdev)*100);
e_min_R2_rem = rem(e_min_R2(1,1),10);
e_max_R2_rem = 10 - rem(e_max_R2(1,1),10);
if e_min_R2 < 0
e_min_R2 = 0;
e_min_R2_rem = 0;
end
if e_max_R2 >= 100
e_max_R2 = 100;
e_max_R2_rem = 0;
end
e_R2 = (e_min_R2-e_min_R2_rem):INT:(e_max_R2+e_max_R2_rem);
end
if N == 5
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2);
end

if N >= 6
f_min_R2 = round((Results_R2_mean_std(6,1) - Results_R2_mean_std(6,2)*stdev)*100);
f_max_R2 = round((Results_R2_mean_std(6,1) + Results_R2_mean_std(6,2)*stdev)*100);
f_min_R2_rem = rem(f_min_R2(1,1),10);
f_max_R2_rem = 10 - rem(f_max_R2(1,1),10);
if f_min_R2 < 0
f_min_R2 = 0;
f_min_R2_rem = 0;
end
if f_max_R2 >= 100
f_max_R2 = 100;
f_max_R2_rem = 0;
end
f_R2 = (f_min_R2-f_min_R2_rem):INT:(f_max_R2+f_max_R2_rem);
end
if N == 6
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2);
end

if N >= 7
g_min_R2 = round((Results_R2_mean_std(7,1) - Results_R2_mean_std(7,2)*stdev)*100);
g_max_R2 = round((Results_R2_mean_std(7,1) + Results_R2_mean_std(7,2)*stdev)*100);
g_min_R2_rem = rem(g_min_R2(1,1),10);
g_max_R2_rem = 10 - rem(g_max_R2(1,1),10);
if g_min_R2 < 0
g_min_R2 = 0;
g_min_R2_rem = 0;
end
if g_max_R2 >= 100
g_max_R2 = 100;
g_max_R2_rem = 0;
end
g_R2 = (g_min_R2-g_min_R2_rem):INT:(g_max_R2+g_max_R2_rem);
end
if N == 7
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2);
end

if N >= 8
h_min_R2 = round((Results_R2_mean_std(8,1) - Results_R2_mean_std(8,2)*stdev)*100);
h_max_R2 = round((Results_R2_mean_std(8,1) + Results_R2_mean_std(8,2)*stdev)*100);
h_min_R2_rem = rem(h_min_R2(1,1),10);
h_max_R2_rem = 10 - rem(h_max_R2(1,1),10);
if h_min_R2 < 0
h_min_R2 = 0;
h_min_R2_rem = 0;
end
if h_max_R2 >= 100
h_max_R2 = 100;
h_max_R2_rem = 0;
end
h_R2 = (h_min_R2-h_min_R2_rem):INT:(h_max_R2+h_max_R2_rem);
end
if N == 8
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2);
end

if N >= 9
l_min_R2 = round((Results_R2_mean_std(9,1) - Results_R2_mean_std(9,2)*stdev)*100);
l_max_R2 = round((Results_R2_mean_std(9,1) + Results_R2_mean_std(9,2)*stdev)*100);
l_min_R2_rem = rem(l_min_R2(1,1),10);
l_max_R2_rem = 10 - rem(l_max_R2(1,1),10);
if l_min_R2 < 0
l_min_R2 = 0;
l_min_R2_rem = 0;
end
if l_max_R2 >= 100
l_max_R2 = 100;
l_max_R2_rem = 0;
end
l_R2 = (l_min_R2-l_min_R2_rem):INT:(l_max_R2+l_max_R2_rem);
end
if N == 9
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2,l_R2);
end

if N >= 10
q_min_R2 = round((Results_R2_mean_std(10,1) - Results_R2_mean_std(10,2)*stdev)*100);
q_max_R2 = round((Results_R2_mean_std(10,1) + Results_R2_mean_std(10,2)*stdev)*100);
q_min_R2_rem = rem(q_min_R2(1,1),10);
q_max_R2_rem = 10 - rem(q_max_R2(1,1),10);
if q_min_R2 < 0
q_min_R2 = 0;
q_min_R2_rem = 0;
end
if q_max_R2 >= 100
q_max_R2 = 100;
q_max_R2_rem = 0;
end
q_R2 = (q_min_R2-q_min_R2_rem):INT:(q_max_R2+q_max_R2_rem);
end
if N == 10
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2,l_R2,q_R2);
end

weights_R2 = weights_allcomb_R2;
number_of_weights_R2=length(weights_allcomb_R2);

D_new = [];
D_new_out = [];
D_conc_out = [];
V_new = [];
V_new_out = [];
V_conc_out = [];
R2_new = [];
R2_new_out = [];
R2_conc_out = [];

end % End if t = 1

% If min and/or max cutoffs are odd numbers after iteration 2 (5% grid spacing) then widen to even numbers so divisible by 2
if t == 3
for g = 1:N
if rem(contrib_min_D(1,g),2) == 1
contrib_min_D(1,g) = contrib_min_D(1,g)-1;
end
end
for g = 1:N
if rem(contrib_max_D(1,g),2) == 1
contrib_max_D(1,g) = contrib_max_D(1,g)+1;
end
end
for g = 1:N
if rem(contrib_min_V(1,g),2) == 1
contrib_min_V(1,g) = contrib_min_V(1,g)-1;
end
end
for g = 1:N
if rem(contrib_max_V(1,g),2) == 1
contrib_max_V(1,g) = contrib_max_V(1,g)+1;
end
end
for g = 1:N
if rem(contrib_min_R2(1,g),2) == 1
contrib_min_R2(1,g) = contrib_min_R2(1,g)-1;
end
end
for g = 1:N
if rem(contrib_max_R2(1,g),2) == 1
contrib_max_R2(1,g) = contrib_max_R2(1,g)+1;
end
end
end

% Combine all possible weight combinations then keep only those that sum to 100%
if t > 1

a_D = contrib_min_D(1,1):INT:contrib_max_D(1,1);
b_D = contrib_min_D(1,2):INT:contrib_max_D(1,2);
if N == 2
weights_allcomb_D = allcomb(a_D,b_D);
end

if N >= 3 
c_D = contrib_min_D(1,3):INT:contrib_max_D(1,3);
end
if N == 3
weights_allcomb_D = allcomb(a_D,b_D,c_D);
end

if N >= 4
d_D = contrib_min_D(1,4):INT:contrib_max_D(1,4);
end
if N == 4
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D);
end

if N >= 5
e_D = contrib_min_D(1,5):INT:contrib_max_D(1,5);
end
if N == 5
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D);
end

if N >= 6
f_D = contrib_min_D(1,6):INT:contrib_max_D(1,6);
end
if N == 6
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D);
end

if N >= 7
g_D = contrib_min_D(1,7):INT:contrib_max_D(1,7);
end
if N == 7
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D);
end

if N >= 8
h_D = contrib_min_D(1,8):INT:contrib_max_D(1,8);
end
if N == 8
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D);
end

if N >= 9
l_D = contrib_min_D(1,9):INT:contrib_max_D(1,9);
end
if N == 9
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D,l_D);
end

if N >= 10
q_D = contrib_min_D(1,10):INT:contrib_max_D(1,10);
end
if N == 10
weights_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D,l_D,q_D);
end

a_V = contrib_min_V(1,1):INT:contrib_max_V(1,1);
b_V = contrib_min_V(1,2):INT:contrib_max_V(1,2);
if N == 2
weights_allcomb_V = allcomb(a_V,b_V);
end

if N >= 3 
c_V = contrib_min_V(1,3):INT:contrib_max_V(1,3);
end
if N == 3
weights_allcomb_V = allcomb(a_V,b_V,c_V);
end

if N >= 4
d_V = contrib_min_V(1,4):INT:contrib_max_V(1,4);
end
if N == 4
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V);
end

if N >= 5
e_V = contrib_min_V(1,5):INT:contrib_max_V(1,5);
end
if N == 5
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V);
end

if N >= 6
f_V = contrib_min_V(1,6):INT:contrib_max_V(1,6);
end
if N == 6
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V);
end

if N >= 7
g_V = contrib_min_V(1,7):INT:contrib_max_V(1,7);
end
if N == 7
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V);
end

if N >= 8
h_V = contrib_min_V(1,8):INT:contrib_max_V(1,8);
end
if N == 8
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V);
end

if N >= 9
l_V = contrib_min_V(1,9):INT:contrib_max_V(1,9);
end
if N == 9
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V,l_V);
end

if N >= 10
q_V = contrib_min_V(1,10):INT:contrib_max_V(1,10);
end
if N == 10
weights_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V,l_V,q_V);
end

a_R2 = contrib_min_R2(1,1):INT:contrib_max_R2(1,1);
b_R2 = contrib_min_R2(1,2):INT:contrib_max_R2(1,2);
if N == 2
weights_allcomb_R2 = allcomb(a_R2,b_R2);
end

if N >= 3 
c_R2 = contrib_min_R2(1,3):INT:contrib_max_R2(1,3);
end
if N == 3
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2);
end

if N >= 4
d_R2 = contrib_min_R2(1,4):INT:contrib_max_R2(1,4);
end
if N == 4
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2);
end

if N >= 5
e_R2 = contrib_min_R2(1,5):INT:contrib_max_R2(1,5);
end
if N == 5
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2);
end

if N >= 6
f_R2 = contrib_min_R2(1,6):INT:contrib_max_R2(1,6);
end
if N == 6
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2);
end

if N >= 7
g_R2 = contrib_min_R2(1,7):INT:contrib_max_R2(1,7);
end
if N == 7
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2);
end

if N >= 8
h_R2 = contrib_min_R2(1,8):INT:contrib_max_R2(1,8);
end
if N == 8
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2);
end

if N >= 9
l_R2 = contrib_min_R2(1,9):INT:contrib_max_R2(1,9);
end
if N == 9
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2,l_R2);
end

if N >= 10
q_R2 = contrib_min_R2(1,10):INT:contrib_max_R2(1,10);
end
if N == 10
weights_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2,l_R2,q_R2);
end

end

% Only keep groups of weights that sum to 100%
C_D = [];
number_of_weights_D = [];
weights_D = [];
C_D=sum(weights_allcomb_D,2);
number_of_weights_D=length(find(C_D==100));
weights_D(1:number_of_weights_D,:) = weights_allcomb_D(find(C_D==100),:);
clear weights_allcomb_D

C_V = [];
number_of_weights_V = [];
weights_V = [];
C_V=sum(weights_allcomb_V,2);
number_of_weights_V=length(find(C_V==100));
weights_V(1:number_of_weights_V,:) = weights_allcomb_V(find(C_V==100),:);
clear weights_allcomb_V

C_R2 = [];
number_of_weights_R2 = [];
weights_R2 = [];
C_R2=sum(weights_allcomb_R2,2);
number_of_weights_R2=length(find(C_R2==100));
weights_R2(1:number_of_weights_R2,:) = weights_allcomb_R2(find(C_R2==100),:);
clear weights_allcomb_R2


if t == 1
weights_10 = weights_R2;
end
if t == 2
weights_05 = weights_R2;
end
if t == 3
weights_02 = weights_R2;
end
if t == 4
weights_01 = weights_R2;
end






%%%%%%%%%%%%%%%%%%%%%%
% ITERATIONS APPLIED %
%%%%%%%%%%%%%%%%%%%%%%

sw_D = [];
sw_V = [];
sw_R2 = [];
D_OPTIMIZE1 = [];
V_OPTIMIZE1 = [];
R2_OPTIMIZE1 = [];

count = 0;
if length(weights_D) > 0
for k = 1:length(weights_D(:,1))
count = count + 1;
cdf_OPTIMIZE1 = sum(cdf_source.*repmat(weights_D(k,:).*0.01,length(x_all),1),2);
D_OPTIMIZE1(count,1) = max([max(cdf_sink - cdf_OPTIMIZE1)],[max(cdf_OPTIMIZE1 - cdf_sink)]);
end % End weights loop
end % End if number of weights for iteration is > 0

count = 0;
if length(weights_V) > 0
for k = 1:length(weights_V(:,1))
count = count + 1;
cdf_OPTIMIZE1 = sum(cdf_source.*repmat(weights_V(k,:).*0.01,length(x_all),1),2);
V_OPTIMIZE1(count,1) = max(cdf_sink - cdf_OPTIMIZE1) + max(cdf_OPTIMIZE1 - cdf_sink);
end % End weights loop
end % End if number of weights for iteration is > 0

count = 0;
if length(weights_R2) > 0
for k = 1:length(weights_R2(:,1))
count = count + 1;
pdp_OPTIMIZE1 = sum(pdp_source.*repmat(weights_R2(k,:).*0.01,length(x),1),2);
R2_OPTIMIZE1(count,1) = ((sum((pdp_sink - mean(pdp_sink)).*(pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1)).*(pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1)))))))*...
	((sum((pdp_sink - mean(pdp_sink)).*(pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1)).*(pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1)))))));
end % End weights loop
end % End if number of weights for iteration is > 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORT AND APPEND RESULTS, SET BOUNDS ON NEXT ITERATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sw_D = size(weights_D);
if sw_D(1,1) > 0
D_conc = [D_OPTIMIZE1,weights_D];
D_conc_sort = sortrows(D_conc, 1); % sorted best fits 
D_conc_size_1 = size(D_conc_out);
D_conc_size_2 = size(D_conc);
D_conc_out(D_conc_size_1(1,1)+1:(D_conc_size_2(1,1)+D_conc_size_1(1,1)),:) = D_conc;
D_conc_out = sortrows(D_conc_out,1);
D_new_size_1 = size(D_new_out);
if sw_D(1,1) < Max_best_fits 
D_new = D_conc_sort;
else
D_new = D_conc_sort(1:Max_best_fits, :);
end
D_new_size_2 = size(D_new);
D_new_out(D_new_size_1(1,1)+1:(D_new_size_2(1,1)+D_new_size_1(1,1)),:) = D_new;
for i = 1:N
contrib_min_D(1,i) = min(D_new(:,i+1));
contrib_max_D(1,i) = max(D_new(:,i+1));
end
end

sw_V = size(weights_V);
if sw_V(1,1) > 0
V_conc = [V_OPTIMIZE1,weights_V];
V_conc_sort = sortrows(V_conc, 1); % sorted best fits 
V_conc_size_1 = size(V_conc_out);
V_conc_size_2 = size(V_conc);
V_conc_out(V_conc_size_1(1,1)+1:(V_conc_size_2(1,1)+V_conc_size_1(1,1)),:) = V_conc;
V_conc_out = sortrows(V_conc_out,1);
V_new_size_1 = size(V_new_out);
if sw_V(1,1) < Max_best_fits 
V_new = V_conc_sort;
else
V_new = V_conc_sort(1:Max_best_fits, :);
end
V_new_size_2 = size(V_new);
V_new_out(V_new_size_1(1,1)+1:(V_new_size_2(1,1)+V_new_size_1(1,1)),:) = V_new;
for i = 1:N
contrib_min_V(1,i) = min(V_new(:,i+1));
contrib_max_V(1,i) = max(V_new(:,i+1));
end
end

sw_R2 = size(weights_R2);
if sw_R2(1,1) > 0
R2_conc = [R2_OPTIMIZE1,weights_R2];
R2_conc_sort = flipud(sortrows(R2_conc, 1)); % sorted best fits 
R2_conc_size_1 = size(R2_conc_out);
R2_conc_size_2 = size(R2_conc);
R2_conc_out(R2_conc_size_1(1,1)+1:(R2_conc_size_2(1,1)+R2_conc_size_1(1,1)),:) = R2_conc;
R2_conc_out = flipud(sortrows(R2_conc_out,1));
R2_new_size_1 = size(R2_new_out);
if sw_R2(1,1) < Max_best_fits 
R2_new = R2_conc_sort;
else
R2_new = R2_conc_sort(1:Max_best_fits, :);
end
R2_new_size_2 = size(R2_new);
R2_new_out(R2_new_size_1(1,1)+1:(R2_new_size_2(1,1)+R2_new_size_1(1,1)),:) = R2_new;
for i = 1:N
contrib_min_R2(1,i) = min(R2_new(:,i+1));
contrib_max_R2(1,i) = max(R2_new(:,i+1));
end
end

contrib_min_D_out(t,:) = contrib_min_D;
contrib_max_D_out(t,:) = contrib_max_D;
contrib_min_V_out(t,:) = contrib_min_V;
contrib_max_V_out(t,:) = contrib_max_V;
contrib_min_R2_out(t,:) = contrib_min_R2;
contrib_max_R2_out(t,:) = contrib_max_R2;

end % End iterations loop

% Calculate mean and stdev. of user-specified Max_best_fits number of best fits and best overall fit
D_mean_OPTIMIZE1 = mean(D_conc_out(1:Max_best_fits,1));
D_std_OPTIMIZE1 = std(D_conc_out(1:Max_best_fits,1));
V_mean_OPTIMIZE1 = mean(V_conc_out(1:Max_best_fits,1));
V_std_OPTIMIZE1 = std(V_conc_out(1:Max_best_fits,1));
R2_mean_OPTIMIZE1 = mean(R2_conc_out(1:Max_best_fits,1));
R2_std_OPTIMIZE1 = std(R2_conc_out(1:Max_best_fits,1));

% Calculate mean and stdev. of best fits and best overall fit. Print results in MATLAB command window
mean_stdev_OPTIMIZE1_D_V_R2 = [[D_mean_OPTIMIZE1; D_std_OPTIMIZE1], [V_mean_OPTIMIZE1; V_std_OPTIMIZE1], [R2_mean_OPTIMIZE1; R2_std_OPTIMIZE1]]
best_OPTIMIZE1_D_V_R2 = [D_conc_out(1,1), V_conc_out(1,1), R2_conc_out(1,1)]
best_OPTIMIZE1_weights = [D_conc_out(1,2:end); V_conc_out(1,2:end); R2_conc_out(1,2:end)].*.01

% Plot all results against mixed sample
%for i = 1:Max_best_fits
cdf_D_best(:,1) = sum(cdf_source.*repmat(D_conc_sort(1,2:end).*0.01,length(x_all),1),2);
%end
figure(F1);
hold on
p2 = plot(x_all,cdf_sink, 'k', 'linewidth', 2);
%for i = 1:Max_best_fits
p1 = plot(x_all,cdf_D_best(:,1), 'color', [0 1 0],'linewidth',1);
%end
hold off
axis([PDP_min PDP_max 0 1])
title('KS test D statistic')
xlabel('Age (Ma)')
ylabel('Cumulative probability');
legend([p1 p2],{'Models','Target'})

%for i = 1:Max_best_fits
cdf_V_best(:,1) = sum(cdf_source.*repmat(V_conc_sort(1,2:end).*0.01,length(x_all),1),2);
%end
figure(F2);
hold on
p4 = plot(x_all,cdf_sink, 'k', 'linewidth', 2);
%for i = 1:Max_best_fits
p3 = plot(x_all,cdf_V_best(:,1), 'color', [0 1 0],'linewidth',1);
%end
hold off
axis([PDP_min PDP_max 0 1])
title('Kuiper test V statistic')
xlabel('Age (Ma)')
ylabel('Cumulative probability');
legend([p3 p4],{'Models','Target'})

%for i = 1:Max_best_fits
pdp_R2_best(:,1) = sum(pdp_source.*repmat(R2_conc_sort(1,2:end).*0.01,length(x),1),2);
%end
figure(F3);
hold on
p6 = plot(x,pdp_sink, 'k', 'linewidth', 2);
%for i = 1:Max_best_fits
p5 = plot(x,pdp_R2_best(:,1), 'color', [0 1 0],'linewidth',1);
%end
hold off
axis([PDP_min PDP_max 0 max(pdp_sink)+max(pdp_sink)*.25])
title('Cross-correlation coefficient')
xlabel('Age (Ma)')
ylabel('Relative probability');
legend([p5 p6],{'Models','Target'})

% Visualize Monte Carlo plots (remove curly brackets to plot)
%{
KS_D = figure;
colours = colormap(jet((N)));
colorbar;
D_conc_out = flipud(D_conc_out);
hold on
for j = 1:length(D_conc_out(:,1))
for i = 1:N
line(D_conc_out(j,i+1),D_conc_out(j,1),'LineStyle','none','Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
	'MarkerSize',6);
end
end
D_conc_out = flipud(D_conc_out);
figure(KS_D)
set(0, 'CurrentFigure', KS_D)
set(KS_D,'Visible','on')
ylabel('KS test D value');
xlabel('Source contribution');
title('Monte Carlo mixture modeling using KS test D value');
set(gca,'Ydir','reverse')
hold off

Kuiper_V = figure;
colours = colormap(jet((N)));
colorbar;
hold on
for j = 1:length(V_conc_out(:,1))
for i = 1:N
V_conc_out = flipud(V_conc_out);
line(V_conc_out(j,i+1),V_conc_out(j,1),'LineStyle','none','Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
	'MarkerSize',6);
V_conc_out = flipud(V_conc_out);
end
end
figure(Kuiper_V)
set(0, 'CurrentFigure', Kuiper_V)
set(Kuiper_V,'Visible','on')
ylabel('Kuiper test V value');
xlabel('Source contribution');
title('Monte Carlo mixture modeling using Kuiper test V value');
set(gca,'Ydir','reverse')
hold off

R2 = figure;
colours = colormap(jet((N)));
colorbar;
hold on
for j = 1:length(R2_conc_out(:,1))
for i = 1:N
R2_conc_out = flipud(R2_conc_out);
line(R2_conc_out(j,i+1),R2_conc_out(j,1),'LineStyle','none','Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
	'MarkerSize',6);
R2_conc_out = flipud(R2_conc_out);
end
end
figure(R2)
set(R2,'Visible','on')
ylabel('Cross-correlation coefficient');
xlabel('Source contribution');
title('Monte Carlo mixture modeling using Cross-correlation coefficient');
hold off
%}

end % End if Optimization routine #1 

%% OPTION TO OPTIMIZE BASED ON BEST FIT MONTE CARLO METHOD UNMIXING RESULTS %%

if Optimize == 2

x = PDP_min:PDP_step:PDP_max;
for i = 1:N+1;
m = data(:,i*2-1);
m = m(isfinite(m(:,1)),:);
s = data(:,i*2);
s = s(isfinite(s(:,1)),:);
f = zeros(length(m),length(x));
for j = 1:length(m);
f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*PDP_step);
end
pdps(:,i) = ((sum(f, 1))/length(m)).';
end
pdp_sink = pdps(:,1); % Mixed sample PDP
pdp_source = pdps(:,2:end); % All source sample PDPs

% Make cumulative distribution functions
binCounts = [];
sumCounts = [];
sumCounts = [];
CDF = [];
for i = 1:N+1
ages(:,i) = (data(:,i*2-1));
end
x_all = sort(nonzeros(reshape(ages,[numel(ages),1])));
x_all(isnan(x_all(:,1)),:) = [];
x_all(x_all>4500) = []; % Remove ages older than Earth
x_all(x_all<0) = []; % Remove future ages
binEdges    =  [-inf ; x_all ; inf];
for i=1:N+1
x1 = ages(:,i);
x1(isnan(x1(:,1)),:) = [];
binCounts(:,i)  =  histc(nonzeros(x1), binEdges, 1);
clear x1
end
for i=1:N+1
sumCounts(:,i)  =  cumsum(binCounts(:,i))./sum(binCounts(:,i));
end
clear binCounts 
for i=1:N+1
CDF(:,i)  =  sumCounts(1:end-1, i);
end
clear sumCounts
cdf_sink = CDF(2:end,1); % Mixed sample CDF
cdf_source = CDF(2:end,2:N+1); % All source sample CDFs
clear CDF

% Get best fits (mean values) from Monte Carlo results. These will be 'initial guesses' in the constrained minimum find algorithm below (fmincon).
weights_mc_D_tmp = sortrows([Passed_D, Passed_weights_D],1);
weights_mc_D = weights_mc_D_tmp(1:Max_best_fits,2:N+1).*100;
weights_mc_V_tmp = sortrows([Passed_V, Passed_weights_V],1);
weights_mc_V = weights_mc_V_tmp(1:Max_best_fits,2:N+1).*100;
weights_mc_R2_tmp = flipud(sortrows([Passed_R2, Passed_weights_R2],1));
weights_mc_R2 = weights_mc_R2_tmp(1:Max_best_fits,2:N+1).*100;

lb = zeros(1,N); % Set lower bound (0%) for possible source contribution 
ub = ones(1,N)*100; % Set upper bound (100%) for possible source contribution 
opts1=  optimset('display','off');

count = 0;

% Apply constrained minimum algorithm. Each metric (D, V, and R2) calls a different function file (e.g., fmincon_D for KS D) that will use initial guesses (best mean values from Monte Carlo results).
for i = 1:Max_best_fits
count = count + 1;

[aa_D, bb_D] = fmincon(@fmincon_D, weights_mc_D(i,:), [], [], [], [], lb, ub, @nonlcon2, opts1);
ba_conc_D(i,:) = [bb_D,aa_D];

[aa_V, bb_V] = fmincon(@fmincon_V, weights_mc_V(i,:), [], [], [], [], lb, ub, @nonlcon2, opts1);
ba_conc_V(i,:) = [bb_V,aa_V];

[aa_R2, bb_R2] = fmincon(@fmincon_R2, weights_mc_R2(i,:), [], [], [], [], lb, ub, @nonlcon2, opts1);
ba_conc_R2(i,:) = [1-bb_R2,aa_R2];

percent_complete = round((count)/Max_best_fits*100)

end

% Concatenate results
ba_sort_D = sortrows(ba_conc_D,1);
ba_sort_V = sortrows(ba_conc_V,1);
ba_sort_R2 = flipud(sortrows(ba_conc_R2,1));

% Sort and trim results 
D_mean_OPTIMIZE2 = mean(ba_sort_D(1:Max_best_fits,1));
D_std_OPTIMIZE2 = std(ba_sort_D(1:Max_best_fits,1));
V_mean_OPTIMIZE2 = mean(ba_sort_V(1:Max_best_fits,1));
V_std_OPTIMIZE2 = std(ba_sort_V(1:Max_best_fits,1));
R2_mean_OPTIMIZE2 = mean(ba_sort_R2(1:Max_best_fits,1));
R2_std_OPTIMIZE2 = std(ba_sort_R2(1:Max_best_fits,1));

% Calculate mean and stdev. of best fits and best overall fit. Print results in MATLAB command window
mean_stdev_OPTIMIZE2_D_V_R2 = [[D_mean_OPTIMIZE2; D_std_OPTIMIZE2], [V_mean_OPTIMIZE2; V_std_OPTIMIZE2], [R2_mean_OPTIMIZE2; R2_std_OPTIMIZE2]]
best_OPTIMIZE2_D_V_R2 = [ba_sort_D(1,1), ba_sort_V(1,1), ba_sort_R2(1,1)]
best_OPTIMIZE2_weights = [ba_sort_D(1,2:end); ba_sort_V(1,2:end); ba_sort_R2(1,2:end)].*.01

% Plot all results
for i = 1:Max_best_fits
cdf_D_best(:,i) = sum(cdf_source.*repmat(ba_sort_D(i,2:end).*0.01,length(x_all),1),2);
end
figure(F1);
hold on
p2 = plot(x_all,cdf_sink, 'k', 'linewidth', 3);
for i = 1:Max_best_fits
p1 = plot(x_all,cdf_D_best(:,i), 'color', [0 1 0],'linewidth',1);
end
hold off
axis([PDP_min PDP_max 0 1])
title('KS test D statistic')
xlabel('Age (Ma)')
ylabel('Cumulative probability');
legend([p1 p2],{'Models','Target'})

for i = 1:Max_best_fits
cdf_V_best(:,i) = sum(cdf_source.*repmat(ba_sort_V(i,2:end).*0.01,length(x_all),1),2);
end
figure(F2);
hold on
p4 = plot(x_all,cdf_sink, 'k', 'linewidth', 3);
for i = 1:Max_best_fits
p3 = plot(x_all,cdf_V_best(:,i), 'color', [0 1 0],'linewidth',1);
end
hold off
axis([PDP_min PDP_max 0 1])
title('Kuiper test V statistic')
xlabel('Age (Ma)')
ylabel('Cumulative probability');
legend([p3 p4],{'Models','Target'})

for i = 1:Max_best_fits
pdp_R2_best(:,i) = sum(pdp_source.*repmat(ba_sort_R2(i,2:end).*0.01,length(x),1),2);
end
figure(F3);
hold on
p6 = plot(x,pdp_sink, 'k', 'linewidth', 3);
for i = 1:Max_best_fits
p5 = plot(x,pdp_R2_best(:,i), 'color', [0 1 0],'linewidth',1);
end
hold off
axis([PDP_min PDP_max 0 max(pdp_sink)+max(pdp_sink)*.25])
title('Cross-correlation coefficient')
xlabel('Age (Ma)')
ylabel('Relative probability');
legend([p5 p6],{'Models','Target'})

end % End if Optimization routine #2
toc