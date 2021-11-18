function varargout = DZmix(varargin)
%DZMIX M-file for DZmix.fig
%      DZMIX, by itself, creates a new DZMIX or raises the existing
%      singleton*.
%
%      H = DZMIX returns the handle to a new DZMIX or the handle to
%      the existing singleton*.
%
%      DZMIX('Property','Value',...) creates a new DZMIX using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to DZmix_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DZMIX('CALLBACK') and DZMIX('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DZMIX.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DZmix

% Last Modified by GUIDE v2.5 25-Apr-2017 10:36:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
	'gui_Singleton',  gui_Singleton, ...
	'gui_OpeningFcn', @DZmix_OpeningFcn, ...
	'gui_OutputFcn',  @DZmix_OutputFcn, ...
	'gui_LayoutFcn',  [], ...
	'gui_Callback',   []);
if nargin && ischar(varargin{1})
	gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
	[varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
	gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DZmix is made visible.
function DZmix_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
imshow('uhlogo_red.png', 'Parent', handles.axes1);
% Choose default command line output for DZmix
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DZmix wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DZmix_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% Step 1 Select File %%
% --- Executes on button press in Browser.
function Browser_Callback(hObject, eventdata, handles)
% hObject    handle to Browser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Open browser window to find and load sample data; data need to be organized into age-uncertainty pairs (ages, uncertainties, ages, uncertianties, ...) with headers (sample names must start with a letter); first 2 columns are the mixed sample, all others are sources
cla(handles.CDFs,'reset');
cla(handles.PDPs,'reset');

%Get user inputs
PDP_min = str2num(get(handles.PDP_min,'String'));
PDP_max = str2num(get(handles.PDP_max,'String'));
PDP_step = str2num(get(handles.PDP_step,'String'));
x_axis_min = PDP_min;
x_axis_max = PDP_max;

global data
global N
global x
global x_all
global cdf_sink
%global cdf_source
global pdp_sink
global pdp_source

[filename, pathname] = uigetfile({'*'},'File Selector');
fullpathname = strcat(pathname, filename);
set(handles.Filepath, 'String', fullpathname); %show path name
[numbers, text1, data_tmp] = xlsread(fullpathname);
data = cell2mat(data_tmp(2:end,:));
[dataR,dataC]=size(data);
N = (dataC/2)-1; % Number of source samples
%% calculate smallest sample size
for j=1:dataC
	small_sample(:,j)=size(data(~isnan(data(:,j)),j),1);
end
small_sample=min(small_sample);
nsamples=(dataC/2);
set(handles.numsamples, 'String', nsamples-1);
set(handles.target_n, 'String', 1);
set(handles.all_n, 'String', nsamples);
set(handles.smallest_text,'String', small_sample);
N = nsamples-1;
set(handles.num_ages,'String', small_sample);

%calculate 1 or 2 sigma
rad_on=get(handles.sigma,'selectedobject');
switch rad_on
	case handles.sigma1
		%
	case handles.sigma2
		for i=1:nsamples
			data(:,2*i)=data(:,2*i)./2;
		end
	otherwise
		set(handles.edit_radioselect,'string','');
end

% Make probability density plots
x = PDP_min:PDP_step:PDP_max;
%select distribution
dist_on=get(handles.comparison,'selectedobject');
switch dist_on
	case handles.density_PDPs
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
			clear f
		end
	case handles.density_KDEs
		pdp_out=zeros(size(x,2),nsamples+1);
		bandwidth = str2num(get(handles.kde_bandwidth,'String'))
		pdp_out(:,1) = x;
		for i = 1:nsamples;
			m = data(:,i*2-1);
			m = m(isfinite(m(:,1)),:);
			test=isfinite(bandwidth)
			if test==1;
				s = ones(length(m)).*bandwidth;
				f = zeros(length(m),length(x));
				for j = 1:length(m);
					f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*PDP_step);
				end
				kdeAi = ((sum(f, 1))/length(m)).';
				clear f
				xmesh1i=x;
				bandwidthi=bandwidth;
			else
				[bandwidthi,kdeAi,xmesh1i,cdfi]=kde(m,length(x),PDP_min,PDP_max);
			end
			pdpi=transpose(interp1(xmesh1i, kdeAi, x));
			bandwidth_out(i,2) = bandwidthi;
			assignin('base','pdpi',pdpi);
			assignin('base','pdp_out',pdp_out);
			pdp_out(:,i+1) = pdpi(:,1);
			pdp_cdfi = transpose(pdpi);
			pdp_normi = pdp_cdfi/sum(pdp_cdfi);
			cumsumi = cumsum(pdp_normi);
			pdp_cdf_out(:,i+1) = (cumsumi);
			pdp_cdf_out(:,1) = x;
		end
		%pdp_out(1:PDP_max, 2:end)=0;
		pdps = pdp_out(:,2:end);
		F = max(pdp_out(:,2:nsamples+1));
		F=max(F);
		
		
	otherwise
		%
end
pdp_sink = pdps(:,1); % Mixed sample PDP
pdp_source = pdps(:,2:end); % All source sample PDPs
assignin('base','pdp_sink',pdp_sink)
assignin('base','pdp_source',pdp_source)
assignin('base','pdps',pdps)

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

%plot PDPs
axes(handles.PDPs);
y_max=max(max(pdps,[],2));
axis(handles.PDPs, [0 PDP_max 0 inf]);
hold on
colours = colormap(jet((nsamples-1)));
for i = 1:nsamples-1;
	plot(x,pdp_source(:,i),'color',colours((i),:),'linewidth',1.5);
	grid on
end
q1=plot(x,pdp_sink, 'k', 'linewidth', 2);
legend([q1],{'Mixed sample'})
switch dist_on
	case handles.density_PDPs
		title('Source and Mixed PDPs')
	case handles.density_KDEs
		title('Source and Mixed KDEs')
end
hold off
xlabel('Age (Ma)')
ylabel('Relative probability');
ax.YAxis.Exponent = -2;



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

%plot CDFs
axes(handles.CDFs);
for i = 2:nsamples;
	datai = data(:,i*2-1);
	datai =datai(isfinite(datai(:,1)),:);
	cdf(i) = cdfplot(datai);
	set(cdf(i),'color',colours((i-1),:),'linewidth',1.5);
	hold on;
	grid on;
	xlabel('');
	ylabel('');
end
q2=cdfplot(data(:,1));
set(q2,'color','k','linewidth',2)
axis([PDP_min, PDP_max, 0, 1]);
legend([q2],{'Mixed sample'})
leg_cdf=legend([q2],{'Mixed sample'});
title('Source and Mixed CDFs');
xlabel('Age (Ma)')
ylabel('Cumulative probability');
colorbar;
hold off

assignin('base','x',x);
assignin('base','pdp_source',pdp_source)

handles.small_sample=small_sample;
handles.N=N;
handles.data=data;
handles.nsamples=nsamples;
handles.text1=text1;
handles.pdps=pdps;
handles.cdf_source=cdf_source;
handles.cdf_sink=cdf_sink;
handles.pdp_source=pdp_source;
handles.pdp_sink=pdp_sink;
handles.x_all=x_all;
handles.x=x;
guidata(hObject,handles);

% --- Executes on button press in Mix.
function Mix_Callback(hObject, eventdata, handles)
% hObject    handle to Mix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=handles.data;
nsamples=handles.nsamples;
text1=handles.text1;
N=handles.N;
small_sample=handles.small_sample;
cdf_source=handles.cdf_source;
cdf_sink=handles.cdf_sink;
pdp_source=handles.pdp_source;
pdp_sink=handles.pdp_sink;
x_all=handles.x_all;
x=handles.x;
for i=1:nsamples
	headers(1,i)=text1(1,2*i-1);
end

cla(handles.PDPs,'reset');

age_num = str2num(get(handles.num_ages,'String'));
trials = str2num(get(handles.run,'String'));
percentnum=str2num(get(handles.percentnum,'String'));

plotnum = percentnum*trials/100;
threshold_D = str2num(get(handles.V,'String'));
threshold_V = str2num(get(handles.V,'String'));
threshold_R2 = str2num(get(handles.R2_PDP,'String'));
PDP_min = str2num(get(handles.PDP_min,'String'));
PDP_max = str2num(get(handles.PDP_max,'String'));
PDP_step = str2num(get(handles.PDP_step,'String'));
x_axis_min = PDP_min;
x_axis_max = PDP_max;
%Licht_approach=get(handles.Raw_ages,'value');
Licht_n=str2num(get(handles.num_ages,'String'));
Licht_N=str2num(get(handles.Licht_N,'String'));

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

% Make probability density plots
x = PDP_min:PDP_step:PDP_max;
%select distribution
dist_on=get(handles.comparison,'selectedobject');
switch dist_on
	case handles.density_PDPs
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
			clear f
		end
	case handles.density_KDEs
		pdp_out=zeros(size(x,2),nsamples+1);
		bandwidth = str2num(get(handles.kde_bandwidth,'String'))
		pdp_out(:,1) = x;
		for i = 1:nsamples;
			m = data(:,i*2-1);
			m = m(isfinite(m(:,1)),:);
			test=isfinite(bandwidth)
			if test==1;
				s = ones(length(m)).*bandwidth;
				f = zeros(length(m),length(x));
				for j = 1:length(m);
					f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*PDP_step);
				end
				kdeAi = ((sum(f, 1))/length(m)).';
				clear f
				xmesh1i=x;
				bandwidthi=bandwidth;
			else
				[bandwidthi,kdeAi,xmesh1i,cdfi]=kde(m,length(x),PDP_min,PDP_max);
			end
			pdpi=transpose(interp1(xmesh1i, kdeAi, x));
			bandwidth_out(i,1) = bandwidthi;
			assignin('base','pdpi',pdpi);
			assignin('base','pdp_out',pdp_out);
			pdp_out(:,i+1) = pdpi(:,1);
			pdp_cdfi = transpose(pdpi);
			pdp_normi = pdp_cdfi/sum(pdp_cdfi);
			cumsumi = cumsum(pdp_normi);
			pdp_cdf_out(:,i+1) = (cumsumi);
			pdp_cdf_out(:,1) = x;
		end
		%pdp_out(1:PDP_max, 2:end)=0;
		pdps = pdp_out(:,2:end);
		F = max(pdp_out(:,2:nsamples+1));
		F=max(F);
		set(handles.kde_bandwidth, 'String',(round(10*mean(bandwidth_out))/10));
	otherwise
		%
end

pdp_sink = pdps(:,1); % Mixed sample PDP
pdp_source = pdps(:,2:end); % All source sample PDPs

%plot PDPs
axes(handles.PDPs);
y_max=max(max(pdps,[],2));
axis(handles.PDPs, [0 PDP_max 0 inf]);
hold on
colours = colormap(jet((nsamples-1)));
for i = 1:nsamples-1;
	plot(x,pdp_source(:,i),'color',colours((i),:),'linewidth',1.5);
	grid on
end
q1=plot(x,pdp_sink, 'k', 'linewidth', 2);
legend([q1],{'Mixed sample'})
switch dist_on
	case handles.density_PDPs
		title('Source and Mixed PDPs')
	case handles.density_KDEs
		title('Source and Mixed KDEs')
end
hold off
xlabel('Age (Ma)')
ylabel('Relative probability');
ax.YAxis.Exponent = -2;

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


%% MONTE CARLO MODEL: RANDOMLY GENERATE SETS OF wghtS TO SCALE SOURCE DISTRIBUTIONS FOR COMPARISON TO MIXED SAMPLE DISTRIBUTION %%
count = 0;
h=waitbar(0,'Calculating','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0);
while count < trials % Run Monte Carlo until reach number of user specified trials
	% Check for Cancel button press
	if getappdata(h,'canceling')
		break
	end
	count = count + 1; % Current trial during Monte Carlo model
	
	% Determine randomly generated set of wghts for each trial
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
	wghts_tmp = diff(sampled_y);
	wghts = wghts_tmp(randperm(length(wghts_tmp))); % Randomly generated set of wghts for each trial
	assignin('base','wghts',wghts);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% IF USER-SPECIFIED Licht_approach = 0 THEN RANDOMLY SAMPLED wghtS WILL SCALE ENTIRE SOURCE DISTRIBUTIONS FOR COMPARISON TO MIXED SAMPLE DISTRIBUTION %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%if Licht_approach == 0
	
	rad_on_source_scaling=get(handles.subsample_ages,'Value');
	if rad_on_source_scaling==0
		
		
		% wght and sum all source distributions
		cdf = sum(cdf_source.*(repmat(wghts,length(x_all),1)),2);
		pdp = sum(pdp_source.*(repmat(wghts,length(x),1)),2);
		
		% Compare randomly wghted sources with mixed sample
		D = max([max(cdf_sink - cdf)],[max(cdf - cdf_sink)]);
		V = max(cdf_sink - cdf) + max(cdf - cdf_sink);
		R2 = ((sum((pdp_sink - mean(pdp_sink)).*(pdp - mean(pdp))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp - mean(pdp)).*(pdp - mean(pdp)))))))*...
			((sum((pdp_sink - mean(pdp_sink)).*(pdp - mean(pdp))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp - mean(pdp)).*(pdp - mean(pdp)))))));
		
		D_std = 0;
		V_std = 0;
		R2_std = 0;
		set(handles.Recent_run, 'String', 'Monte Carlo unmixing model (scale source distributions)');
		
		%end % End if Licht_approach = 0
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% IF USER-SPECIFIED Licht_approach = 1 THEN Licht_N SOURCE DISTRIBUTIONS ARE CREATED BY RANDOMLY SAMPLING SOURCE AGES BASED ON trial wght FOR A TOTAL OF Licht_n AGES %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		%if Licht_approach == 1
		
	else
		
		% First determine how many ages correspond to each source based on rounding the randomly determined trial wght
		num_Ages = round(wghts.*Licht_n);
		
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
			break
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
			assignin('base','s',s)
			assignin('base','m',m)
			assignin('base','x',x)
			pdp_rand_ages_all = ((sum(bsxfun(@times,1./(s.*sqrt(2*pi)),exp(bsxfun(@rdivide, -(bsxfun(@minus, x, m).^2), 2*s.^2))), 1)/length(m)).').*PDP_step;
			pdp_source(:,p) = pdp_rand_ages_all;
		end
		pdp = mean(pdp_source,2);
		
		% Compare each randomly generated age distribution to the mixed sample distribution for each trial wght
		for p = 1:Licht_N
			D_tmp(:,p) = max([max(cdf_sink - cdf_source(:,p))],[max(cdf_source(:,p) - cdf_sink)]);
			V_tmp(:,p) = max(cdf_sink - cdf_source(:,p)) + max(cdf_source(:,p) - cdf_sink);
			R2_tmp(:,p) = ((sum((pdp_sink - mean(pdp_sink)).*(pdp_source(:,p) - mean(pdp_source(:,p)))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp_source(:,p) - ...
				mean(pdp_source(:,p))).*(pdp_source(:,p) - mean(pdp_source(:,p))))))))*((sum((pdp_sink - mean(pdp_sink)).*(pdp_source(:,p) - mean(pdp_source(:,p)))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*...
				(pdp_sink - mean(pdp_sink))))*(sum((pdp_source(:,p) - mean(pdp_source(:,p))).*(pdp_source(:,p) - mean(pdp_source(:,p))))))));
		end
		
		% Calculate mean and stdev for each group of comparisons (total of Licht_N) within each trial; will be filtered below based on mean value
		D = mean(D_tmp);
		V = mean(V_tmp);
		R2 = mean(R2_tmp);
		
		D_std = std(D_tmp);
		V_std = std(V_tmp);
		R2_std = std(R2_tmp);
		set(handles.Recent_run, 'String', 'Monte Carlo unmixing model (subsample source ages)');
		
		%end % End if Licht_approach = 1
		
	end % End radio button
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FILTER AND SORT MODEL RESULTS BASED ON USER SPECIFIED THRESHOLDS, DETERMINE BEST FIT SOURCE CONTRIBUTIONS, PLOT RESULTS %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Keep results if they pass user specified thresholds
	if D <= threshold_D
		fit_D = fit_D + 1;
		Passed_cdf_D(:,fit_D) = cdf;
		Passed_wghts_D(fit_D, :) = wghts;
		Passed_D(fit_D, 1) = D;
		Passed_x_all_D(:,fit_D) = x_all;
		Passed_D_std(fit_D, 1) = D_std;
	end
	
	if V <= threshold_V
		fit_V = fit_V + 1;
		Passed_cdf_V(:,fit_V) = cdf;
		Passed_wghts_V(fit_V, :) = wghts;
		Passed_V(fit_V, 1) = V;
		Passed_x_all_V(:,fit_V) = x_all;
		Passed_V_std(fit_V, 1) = V_std;
	end
	
	if R2 >= threshold_R2
		fit_R2 = fit_R2 + 1;
		Passed_pdp_R2(:,fit_R2) = pdp;
		Passed_wghts_R2(fit_R2, :) = wghts;
		Passed_R2(fit_R2, 1) = R2;
		Passed_R2_std(fit_R2,:) = R2_std;
	end
	
	waitbar(count/trials);
end % End once reached number of user specified trials
delete (h)

pdp_R2_export=Passed_pdp_R2(:,1:plotnum);
cdf_D_export=Passed_cdf_D(:,1:plotnum);
cdf_V_export=Passed_cdf_V(:,1:plotnum);

% Sort and trim results to retain 'plotnum' number of best model fits
if length(Passed_D) >= plotnum
	Results_D = sortrows([Passed_D, Passed_D_std, Passed_wghts_D, transpose(Passed_cdf_D), transpose(Passed_x_all_D)], 1);
	Results_D = Results_D(1:plotnum,:);
	Results_D_mean_std = transpose([mean(Results_D(1:plotnum,3:N+2));std(Results_D(1:plotnum,3:N+2))]);
end

if length(Passed_V) >= plotnum
	Results_V = sortrows([Passed_V, Passed_V_std, Passed_wghts_V, transpose(Passed_cdf_V), transpose(Passed_x_all_V)], 1);
	Results_V = Results_V(1:plotnum,:);
	Results_V_mean_std = transpose([mean(Results_V(1:plotnum,3:N+2));std(Results_V(1:plotnum,3:N+2))]);
end

if length(Passed_R2) >= plotnum
	Results_R2 = flipud(sortrows([Passed_R2, Passed_R2_std, Passed_wghts_R2, transpose(Passed_pdp_R2)], 1));
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

assignin('base','D_mean',D_mean);

%if Licht_approach == 0
if rad_on_source_scaling==0;
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
	%end
end % End switch

% Calculate mean and stdev. of best fits and best overall fit. Print results in MATLAB command window
mean_stdev_MONTECARLO_D_V_R2 = [[D_mean; D_std], [V_mean; V_std], [R2_mean; R2_std]];
best_MONTECARLO_D_V_R2 = [min(Passed_D), min(Passed_V), max(Passed_R2)];
best_MONTECARLO_wghts = [Results_D(1,3:N+2); Results_V(1,3:N+2); Results_R2(1,3:N+2)];

set(handles.D_std_dev2, 'String', 'Mean D value:');
set(handles.V_std_dev2, 'String', 'Mean V value:');
set(handles.R2_std_dev2, 'String', 'Mean Cross-correlation coefficient:');
set(handles.plus_D, 'String', '+/-');
set(handles.plus_V, 'String', '+/-');
set(handles.plus_R2, 'String', '+/-');
set(handles.mean_R2, 'String', round(R2_mean*1000)/1000);
set(handles.mean_V,'String',round(V_mean*1000)/1000);
set(handles.mean_D,'String',round(D_mean*1000)/1000);
set(handles.R2_std, 'String', round(R2_std*1000)/1000);
set(handles.V_std,'String',round(V_std*1000)/1000);
set(handles.D_std,'String',round(D_std*1000)/1000);

% Plot 'plotnum' number of best model fits for KS test D value
if length(Results_D(:,1)) >= plotnum
	axes(handles.KS_plot);
	for i = 1:plotnum
		p1=plot(Results_D(i,N+3+length(Passed_cdf_D(:,1)):end),Results_D(i,N+3:N+2+length(Passed_cdf_D(:,1))),...
			'color', [0 .68 .3], 'linewidth', 1);
		hold on
	end
	cdf_sink_ages = nonzeros(data(:,1));
	cdf_sink_ages(isnan(cdf_sink_ages(:,1)),:) = [];
	p2=cdfplot(cdf_sink_ages(:,1));
	set(p2,'color','k','linewidth',2);
	hold off
	axis([PDP_min PDP_max 0 1])
	title('KS test D Statistic')
	xlabel('Age (Ma)')
	ylabel('Cumulative probability');
	legend([p1 p2],{'Models','Mixed sample'})
end

% Plot 'plotnum' number of best model fits for Kuiper test V value
if length(Results_V(:,1)) >= plotnum
	axes(handles.Kuiper_plot)
	for i = 1:plotnum
		p3=plot(Results_V(i,N+3+length(Passed_cdf_V(:,1)):end),Results_V(i,N+3:N+2+length(Passed_cdf_V(:,1))),...
			'color', [1 .1 0], 'linewidth', 1);
		hold on
	end
	cdf_sink_ages = nonzeros(data(:,1));
	cdf_sink_ages(isnan(cdf_sink_ages(:,1)),:) = [];
	p4=cdfplot(cdf_sink_ages(:,1));
	set(p4,'color','k','linewidth',2)
	hold off
	axis([PDP_min PDP_max 0 1])
	title('Kuier test V Statistic')
	xlabel('Age (Ma)')
	ylabel('Cumulative probability');
	legend([p3 p4],{'Models','Mixed sample'})
end

% Plot 'plotnum' number of best model fits for cross correlation (R2) of PDPs
if length(Results_R2) >= plotnum
	axes(handles.R2_plot)
	for i = 1:plotnum
		p5=plot(x,Results_R2(i,N+3:end), 'color', [.04 .31 1], 'linewidth', 1);
		hold on
	end
	p6=plot(x,pdp_sink, 'k', 'linewidth', 2);
	hold off
	axis([PDP_min PDP_max 0 max(pdp_sink)+max(pdp_sink)*.25])
	title('Cross-correlation coefficient')
	xlabel('Age (Ma)')
	ylabel('Relative probability');
	legend([p5 p6],{'Models','Mixed sample'})
end

%update handles for results export before the results are trimmed in
%the next step
handles.Results_D=Results_D;
handles.Results_V=Results_V;
handles.Results_R2=Results_R2;

%output mean and standard deviations of model contributions to tables in
%GUI interface
colnames={'Sample Names', '   Relative Contribution', '  Standard Deviation'};
rows=transpose(headers);
rows=rows(2:end,1);
if length(Passed_D) >= plotnum
	Results_D = sortrows([Passed_D, Passed_wghts_D, transpose(Passed_cdf_D)], 1);
	Results_D = Results_D(1:plotnum,:);
	Results_D_mean_std = transpose([mean(Results_D(1:plotnum,2:N+1));std(Results_D(1:plotnum,2:N+1))]);
	assignin('base','Results_D_mean_std',Results_D_mean_std)
	data_D=num2cell(Results_D_mean_std);
	data_D_table=horzcat(rows,data_D);
	set(handles.D_table,'data',data_D_table, 'ColumnName', colnames);
end
if length(Passed_V) >= plotnum
	Results_V = sortrows([Passed_V, Passed_wghts_V, transpose(Passed_cdf_V)], 1);
	Results_V = Results_V(1:plotnum,:);
	Results_V_mean_std = transpose([mean(Results_V(1:plotnum,2:N+1));std(Results_V(1:plotnum,2:N+1))]);
	data_V=num2cell(Results_V_mean_std);
	data_V_table=horzcat(rows,data_V);
	set(handles.V_table,'data',data_V_table, 'ColumnName', colnames);
end
if length(Passed_R2) >= plotnum
	Results_R2 = flipud(sortrows([Passed_R2, Passed_wghts_R2, transpose(Passed_pdp_R2)], 1));
	Results_R2 = Results_R2(1:plotnum,:);
	Results_R2_mean_std = transpose([mean(Results_R2(1:plotnum,2:N+1));std(Results_R2(1:plotnum,2:N+1))]);
	data_R2=num2cell(Results_R2_mean_std);
	data_R2_table=horzcat(rows,data_R2);
	set(handles.R2_table,'data',data_R2_table, 'ColumnName', colnames)
end

%Concatenate results for export
Results_export= cell(14+3*N,3);
Results_export(1,1)={'Results of Monte Carlo unmixing model'};
%if Licht_approach==0
if rad_on_source_scaling==0;
	
	
	Results_export(2,1)={'wghtings applied to PDPs and CDFs'};
	Results_export(3,:)={'number of grains in each trial=' 'N/A' ''};
else
	
	age_num = str2num(get(handles.num_ages,'String'));
	Results_export(2,1)={'wghtings applied to raw ages'};
	Results_export(3,:)={'number of grains in each trial=' age_num ''};
end % End switch

Results_export(4,:)={'trials=' trials ''};
Results_export(5,:)={'percent accepted trials=' percentnum ''};
Results_export(6,:)={'Cross-correlation cutoff=' threshold_R2 ''};
Results_export(7,:)={'Kuiper V value cutoff=' threshold_V ''};
Results_export(8,:)={'KS D value cutoff=' threshold_D ''};
Results_export(10,:)={'Cross-correlation' '' ''};
Results_export(11,:)=colnames;
for i=1:N
	Results_export(i+11,:)=data_R2_table(i,:);
end
Results_export(13+N,1)={'Kuiper V value'};
Results_export(14+N,:)=colnames;
for i=1:N
	Results_export(i+14+N,:)=data_V_table(i,:);
end
Results_export(16+2*N,1)={'KS D value'};
Results_export(17+2*N,:)=colnames;
for i=1:N
	Results_export(i+17+2*N,:)=data_D_table(i,:);
end

Results_export(10,2) = strcat(get(handles.mean_R2,'string'),{' '},get(handles.plus_R2,'string'),{' '},get(handles.R2_std,'string'));
Results_export(13+N,2) = strcat(get(handles.mean_D,'string'),{' '},get(handles.plus_D,'string'),{' '},get(handles.D_std,'string'));
Results_export(16+2*N,2) = strcat(get(handles.mean_V,'string'),{' '},get(handles.plus_V,'string'),{' '},get(handles.V_std,'string'));

optimization=0;
%update handles
handles.cdf_V_export=cdf_V_export;
handles.cdf_D_export=cdf_D_export;
handles.pdp_R2_export=pdp_R2_export;
handles.cdf_source=cdf_source;
handles.cdf_sink=cdf_sink;
handles.data_D=data_D;
handles.data_V=data_V;
handles.data_R2=data_R2;
handles.headers=headers;
handles.Results_V_mean_std=Results_V_mean_std;
handles.Results_D_mean_std=Results_D_mean_std;
handles.Results_R2_mean_std=Results_R2_mean_std;
handles.x_all=x_all;
handles.x=x;
handles.optimization=optimization;
handles.pdp_source=pdp_source;
handles.pdp_sink=pdp_sink;
handles.Results_export=Results_export;
handles.plotnum=plotnum;
handles.Passed_cdf_D=Passed_cdf_D;
handles.N=N;
handles.Passed_D=Passed_D;
handles.Passed_V=Passed_V;
handles.Passed_R2=Passed_R2;
handles.Passed_wghts_D=Passed_wghts_D;
handles.Passed_wghts_V=Passed_wghts_V;
handles.Passed_wghts_R2=Passed_wghts_R2;
guidata(hObject,handles);



%% --- Executes on button press in Optimize.
function Optimize_Callback(hObject, eventdata, handles)
% hObject    handle to Optimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global pdp_source
%Get user inputs
Max_best_fits = str2num(get(handles.Max_best_fits,'String'));
percentnum=str2num(get(handles.percentnum,'String'));
%percentnum=2.5 %hardwire in top percent accepted
trials = str2num(get(handles.run,'String'));
Licht_N=trials;
plotnum = percentnum*trials/100;
threshold_D = str2num(get(handles.V,'String'));
threshold_V = str2num(get(handles.V,'String'));
threshold_R2 = str2num(get(handles.R2_PDP,'String'));
PDP_min = str2num(get(handles.PDP_min,'String'));
PDP_max = str2num(get(handles.PDP_max,'String'));
PDP_step = str2num(get(handles.PDP_step,'String'));
x_axis_min = PDP_min;
x_axis_max = PDP_max;
PDP_min = str2num(get(handles.PDP_min,'String'));
PDP_max = str2num(get(handles.PDP_max,'String'));
PDP_step = str2num(get(handles.PDP_step,'String'));
x_axis_min = PDP_min;
x_axis_max = PDP_max;
optimization=0;
%Licht_approach=get(handles.Raw_ages,'value');
Licht_n=str2num(get(handles.num_ages,'String'));
global cdf_source
%Get handles
cdf_source=handles.cdf_source;
cdf_sink=handles.cdf_sink;
headers=handles.headers;
N=handles.N;
pdp_source=handles.pdp_source;
pdp_sink=handles.pdp_sink;
Passed_D=handles.Passed_D;
Passed_V=handles.Passed_V;
Passed_R2=handles.Passed_R2;
Passed_wghts_D=handles.Passed_wghts_D;
Passed_wghts_V=handles.Passed_wghts_V;
Passed_wghts_R2=handles.Passed_wghts_R2;
Results_D_mean_std=handles.Results_D_mean_std;
Results_V_mean_std=handles.Results_V_mean_std;
Results_R2_mean_std=handles.Results_R2_mean_std;
x=handles.x;
data=handles.data;

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
cdf_source=[];
cdf_source = CDF(2:end,2:N+1); % All source sample CDFs
clear CDF

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
	clear f
end
pdp_sink = pdps(:,1); % Mixed sample PDP
pdp_source=[];
pdp_source = pdps(:,2:end); % All source sample PDPs

%% OPTIMIZATION OPTION #1, CONSTRAINTS SET BY MEAN AND STANDARD DEVIATION OF BEST FIT MONTE CARLO MODEL RESULTS %%
opt=get(handles.opt,'selectedobject');
switch opt
	case handles.Sundell_optimize
		optimization=1;
		handles.optimization=optimization;
		set(handles.Recent_run, 'String', 'Iterative optimization unmixing model');
		stdev = 1;
		iteration_int = [10, 5, 2, 1]; % Grid spacing (in percent). First find best fits at 10% spacing, then 5%, then 2%, then 1%.
		iter=0;
		w=waitbar(0,'Iterations','CreateCancelBtn',...
			'setappdata(gcbf,''canceling'',1)');
		setappdata(w,'canceling',0);
		for t = 1:length(iteration_int)
			iter=iter+1;
			% Check for Cancel button press
			if getappdata(w,'canceling')
				break
			end
			INT = iteration_int(1,t); % Set grid spacing
			
			nodes = 0; % Reset number of different wght possibilities
			
			if t == 1
				
				% Set initial grid spacing based on Monte Carlo mean and stdev. and force to be divisibl by 10. For example, if mean wght is 22+/-5 (range of 17 to 27), grid for that source would be set at 10, 20, 30 for iteration 1.
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
					wghts_allcomb_D = allcomb(a_D,b_D);
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
					wghts_allcomb_D = allcomb(a_D,b_D,c_D);
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
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D);
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
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D);
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
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D);
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
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D);
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
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D);
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
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D,l_D);
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
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D,l_D,q_D);
				end
				
				wghts_D = wghts_allcomb_D;
				number_of_wghts_D=length(wghts_allcomb_D);
				
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
					wghts_allcomb_V = allcomb(a_V,b_V);
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
					wghts_allcomb_V = allcomb(a_V,b_V,c_V);
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
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V);
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
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V);
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
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V);
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
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V);
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
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V);
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
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V,l_V);
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
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V,l_V,q_V);
				end
				
				wghts_V = wghts_allcomb_V;
				number_of_wghts_V=length(wghts_allcomb_V);
				
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
					wghts_allcomb_R2 = allcomb(a_R2,b_R2);
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
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2);
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
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2);
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
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2);
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
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2);
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
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2);
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
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2);
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
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2,l_R2);
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
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2,l_R2,q_R2);
				end
				
				wghts_R2 = wghts_allcomb_R2;
				number_of_wghts_R2=length(wghts_allcomb_R2);
				
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
			
			% Combine all possible wght combinations then keep only those that sum to 100%
			if t > 1
				
				a_D = contrib_min_D(1,1):INT:contrib_max_D(1,1);
				b_D = contrib_min_D(1,2):INT:contrib_max_D(1,2);
				if N == 2
					wghts_allcomb_D = allcomb(a_D,b_D);
				end
				
				if N >= 3
					c_D = contrib_min_D(1,3):INT:contrib_max_D(1,3);
				end
				if N == 3
					wghts_allcomb_D = allcomb(a_D,b_D,c_D);
				end
				
				if N >= 4
					d_D = contrib_min_D(1,4):INT:contrib_max_D(1,4);
				end
				if N == 4
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D);
				end
				
				if N >= 5
					e_D = contrib_min_D(1,5):INT:contrib_max_D(1,5);
				end
				if N == 5
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D);
				end
				
				if N >= 6
					f_D = contrib_min_D(1,6):INT:contrib_max_D(1,6);
				end
				if N == 6
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D);
				end
				
				if N >= 7
					g_D = contrib_min_D(1,7):INT:contrib_max_D(1,7);
				end
				if N == 7
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D);
				end
				
				if N >= 8
					h_D = contrib_min_D(1,8):INT:contrib_max_D(1,8);
				end
				if N == 8
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D);
				end
				
				if N >= 9
					l_D = contrib_min_D(1,9):INT:contrib_max_D(1,9);
				end
				if N == 9
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D,l_D);
				end
				
				if N >= 10
					q_D = contrib_min_D(1,10):INT:contrib_max_D(1,10);
				end
				if N == 10
					wghts_allcomb_D = allcomb(a_D,b_D,c_D,d_D,e_D,f_D,g_D,h_D,l_D,q_D);
				end
				
				a_V = contrib_min_V(1,1):INT:contrib_max_V(1,1);
				b_V = contrib_min_V(1,2):INT:contrib_max_V(1,2);
				if N == 2
					wghts_allcomb_V = allcomb(a_V,b_V);
				end
				
				if N >= 3
					c_V = contrib_min_V(1,3):INT:contrib_max_V(1,3);
				end
				if N == 3
					wghts_allcomb_V = allcomb(a_V,b_V,c_V);
				end
				
				if N >= 4
					d_V = contrib_min_V(1,4):INT:contrib_max_V(1,4);
				end
				if N == 4
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V);
				end
				
				if N >= 5
					e_V = contrib_min_V(1,5):INT:contrib_max_V(1,5);
				end
				if N == 5
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V);
				end
				
				if N >= 6
					f_V = contrib_min_V(1,6):INT:contrib_max_V(1,6);
				end
				if N == 6
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V);
				end
				
				if N >= 7
					g_V = contrib_min_V(1,7):INT:contrib_max_V(1,7);
				end
				if N == 7
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V);
				end
				
				if N >= 8
					h_V = contrib_min_V(1,8):INT:contrib_max_V(1,8);
				end
				if N == 8
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V);
				end
				
				if N >= 9
					l_V = contrib_min_V(1,9):INT:contrib_max_V(1,9);
				end
				if N == 9
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V,l_V);
				end
				
				if N >= 10
					q_V = contrib_min_V(1,10):INT:contrib_max_V(1,10);
				end
				if N == 10
					wghts_allcomb_V = allcomb(a_V,b_V,c_V,d_V,e_V,f_V,g_V,h_V,l_V,q_V);
				end
				
				a_R2 = contrib_min_R2(1,1):INT:contrib_max_R2(1,1);
				b_R2 = contrib_min_R2(1,2):INT:contrib_max_R2(1,2);
				if N == 2
					wghts_allcomb_R2 = allcomb(a_R2,b_R2);
				end
				
				if N >= 3
					c_R2 = contrib_min_R2(1,3):INT:contrib_max_R2(1,3);
				end
				if N == 3
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2);
				end
				
				if N >= 4
					d_R2 = contrib_min_R2(1,4):INT:contrib_max_R2(1,4);
				end
				if N == 4
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2);
				end
				
				if N >= 5
					e_R2 = contrib_min_R2(1,5):INT:contrib_max_R2(1,5);
				end
				if N == 5
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2);
				end
				
				if N >= 6
					f_R2 = contrib_min_R2(1,6):INT:contrib_max_R2(1,6);
				end
				if N == 6
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2);
				end
				
				if N >= 7
					g_R2 = contrib_min_R2(1,7):INT:contrib_max_R2(1,7);
				end
				if N == 7
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2);
				end
				
				if N >= 8
					h_R2 = contrib_min_R2(1,8):INT:contrib_max_R2(1,8);
				end
				if N == 8
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2);
				end
				
				if N >= 9
					l_R2 = contrib_min_R2(1,9):INT:contrib_max_R2(1,9);
				end
				if N == 9
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2,l_R2);
				end
				
				if N >= 10
					q_R2 = contrib_min_R2(1,10):INT:contrib_max_R2(1,10);
				end
				if N == 10
					wghts_allcomb_R2 = allcomb(a_R2,b_R2,c_R2,d_R2,e_R2,f_R2,g_R2,h_R2,l_R2,q_R2);
				end
				
			end
			
			% Only keep groups of wghts that sum to 100%
			C_D = [];
			number_of_wghts_D = [];
			wghts_D = [];
			C_D=sum(wghts_allcomb_D,2);
			number_of_wghts_D=length(find(C_D==100));
			wghts_D(1:number_of_wghts_D,:) = wghts_allcomb_D(find(C_D==100),:);
			clear wghts_allcomb_D
			
			C_V = [];
			number_of_wghts_V = [];
			wghts_V = [];
			C_V=sum(wghts_allcomb_V,2);
			number_of_wghts_V=length(find(C_V==100));
			wghts_V(1:number_of_wghts_V,:) = wghts_allcomb_V(find(C_V==100),:);
			clear wghts_allcomb_V
			
			C_R2 = [];
			number_of_wghts_R2 = [];
			wghts_R2 = [];
			C_R2=sum(wghts_allcomb_R2,2);
			number_of_wghts_R2=length(find(C_R2==100));
			wghts_R2(1:number_of_wghts_R2,:) = wghts_allcomb_R2(find(C_R2==100),:);
			clear wghts_allcomb_R2
			
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
			h=waitbar(0,'Calculation (1/3: KS)','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
			setappdata(h,'canceling',0);
			movegui(h,[600,200])
			
			assignin('base','cdf_source',cdf_source)
			assignin('base','wghts_D',wghts_D)
			assignin('base','x_all',x_all)
			
			if length(wghts_D) > 0
				for k = 1:length(wghts_D(:,1))
					% Check for Cancel button press
					if getappdata(h,'canceling')
						break
					end
					count = count + 1;
					cdf_OPTIMIZE1 = sum(cdf_source.*repmat(wghts_D(k,:).*0.01,length(x_all),1),2);
					D_OPTIMIZE1(count,1) = max([max(cdf_sink - cdf_OPTIMIZE1)],[max(cdf_OPTIMIZE1 - cdf_sink)]);
					waitbar(count/length(wghts_D(:,1)),h)
				end % End wghts loop
				delete (h)
			end % End if number of wghts for iteration is > 0
			
			count = 0;
			h=waitbar(0,'Calculation (2/3: Kuiper)','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
			setappdata(h,'canceling',0);
			movegui(h,[600,200])
			if length(wghts_V) > 0
				for k = 1:length(wghts_V(:,1))
					% Check for Cancel button press
					if getappdata(h,'canceling')
						break
					end
					count = count + 1;
					cdf_OPTIMIZE1 = sum(cdf_source.*repmat(wghts_V(k,:).*0.01,length(x_all),1),2);
					V_OPTIMIZE1(count,1) = max(cdf_sink - cdf_OPTIMIZE1) + max(cdf_OPTIMIZE1 - cdf_sink);
					waitbar(count/length(wghts_V(:,1)),h)
				end % End wghts loop
				delete (h)
			end % End if number of wghts for iteration is > 0
			
			count = 0;
			h=waitbar(0,'Calculation (3/3: Cross-correlation coefficient)','CreateCancelBtn',...
				'setappdata(gcbf,''canceling'',1)');
			setappdata(h,'canceling',0);
			movegui(h,[600,200])
			if length(wghts_R2) > 0
				for k = 1:length(wghts_R2(:,1))
					% Check for Cancel button press
					if getappdata(h,'canceling')
						break
					end
					count = count + 1;
					pdp_OPTIMIZE1 = sum(pdp_source.*repmat(wghts_R2(k,:).*0.01,length(x),1),2);
					R2_OPTIMIZE1(count,1) = ((sum((pdp_sink - mean(pdp_sink)).*(pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1)).*(pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1)))))))*...
						((sum((pdp_sink - mean(pdp_sink)).*(pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1)).*(pdp_OPTIMIZE1 - mean(pdp_OPTIMIZE1)))))));
					waitbar(count/length(wghts_R2(:,1)),h)
				end % End wghts loop
				delete (h)
			end % End if number of wghts for iteration is > 0
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% SORT AND APPEND RESULTS, SET BOUNDS ON NEXT ITERATION %
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			sw_D = size(wghts_D);
			if sw_D(1,1) > 0
				D_conc = [D_OPTIMIZE1,wghts_D];
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
			
			sw_V = size(wghts_V);
			if sw_V(1,1) > 0
				V_conc = [V_OPTIMIZE1,wghts_V];
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
			
			sw_R2 = size(wghts_R2);
			if sw_R2(1,1) > 0
				R2_conc = [R2_OPTIMIZE1,wghts_R2];
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
			waitbar(iter/4,w)
		end % End iterations loop
		delete(w)
		
		D_top_percent=round(length(D_conc_out)*percentnum/100);
		V_top_percent=round(length(V_conc_out)*percentnum/100);
		R2_top_percent=round(length(R2_conc_out)*percentnum/100);
		
		% Calculate mean and stdev. of user-specified Max_best_fits number of best fits and best overall fit
		D_best_OPTIMIZE1 = D_conc_out(1,:);
		V_best_OPTIMIZE1 = V_conc_out(1,:);
		R2_best_OPTIMIZE1 = R2_conc_out(1,:);
		
		%{
Calculate mean and stdev. of best fits and best overall fit. Print results in MATLAB command window
mean_stdev_OPTIMIZE1_D_V_R2 = [[D_mean_OPTIMIZE1; D_std_OPTIMIZE1], [V_mean_OPTIMIZE1; V_std_OPTIMIZE1], [R2_mean_OPTIMIZE1; R2_std_OPTIMIZE1]]
best_OPTIMIZE1_D_V_R2 = [D_conc_out(1,1), V_conc_out(1,1), R2_conc_out(1,1)]
		%}
		
		best_OPTIMIZE1_wghts = [D_conc_out(1,2:end); V_conc_out(1,2:end); R2_conc_out(1,2:end)].*.01;
		
		set(handles.mean_R2, 'String', round(R2_best_OPTIMIZE1*1000)/1000)
		set(handles.mean_V,'String',round(V_best_OPTIMIZE1*1000)/1000)
		set(handles.mean_D,'String',round(D_best_OPTIMIZE1*1000)/1000)
		set(handles.R2_std, 'String', 'N/A')
		set(handles.V_std,'String','N/A')
		set(handles.D_std,'String','N/A')
		
		%output average and standard deviations of model contributions to tables in
		%GUI interface
		colnames={'Samples', 'Relative Percent Contribution', 'Standard Deviation'};
		rows=transpose(headers);
		rows=rows(2:end,1);
		
		if length(D_conc_out) >= Max_best_fits
			Results_D = D_conc_out(1,:);
			assignin('base','Results_D',Results_D);
			Results_D_mean_std = transpose([Results_D(1,2:N+1);zeros(1,N)]);
			data_D=num2cell(Results_D_mean_std);
			data_D_table=horzcat(rows,data_D);
			set(handles.D_table,'data',data_D_table, 'ColumnName', colnames);
		end
		if length(V_conc_out) >= Max_best_fits
			Results_V = V_conc_out(1,:);
			Results_V_mean_std = transpose([Results_V(1,2:N+1);zeros(1,N)]);
			data_V=num2cell(Results_V_mean_std);
			data_V_table=horzcat(rows,data_V);
			set(handles.V_table,'data',data_V_table, 'ColumnName', colnames);
		end
		if length(R2_conc_out) >= Max_best_fits
			Results_R2 = R2_conc_out(1,:);
			assignin('base','Results_R2',Results_R2);
			Results_R2_mean_std = transpose([Results_R2(1,2:N+1);zeros(1,N)]);
			data_R2=num2cell(Results_R2_mean_std);
			data_R2_table=horzcat(rows,data_R2);
			set(handles.R2_table,'data',data_R2_table, 'ColumnName', colnames);
		end
		
		%Concatenate results for export
		Results_export= cell(14+3*N,3);
		Results_export(1,1)={'Results of Iterative Optimization unmixing model'};
		%if Licht_approach==0
		rad_on_source_scaling=get(handles.subsample_ages,'Value');
		if rad_on_source_scaling==0
			Results_export(2,1)={'weightings applied to PDPs and CDFs'};
			Results_export(3,:)={'number of grains in each trial=' 'N/A' ''};
		else
			
			age_num = str2num(get(handles.num_ages,'String'));
			Results_export(2,1)={'weightings applied to raw ages'};
			Results_export(3,:)={'number of grains in each trial=' age_num ''};
		end % End switch
		Results_export(4,:)={'trials=' trials ''};
		Results_export(5,:)={'percent accepted trials=' percentnum ''};
		Results_export(6,:)={'Cross-correlation cutoff=' threshold_R2 ''};
		Results_export(7,:)={'Kuiper V value cutoff=' threshold_V ''};
		Results_export(8,:)={'KS D value cutoff=' threshold_D ''};
		Results_export(10,:)={'Cross-correlation' '' ''};
		Results_export(11,:)=colnames;
		for i=1:N
			Results_export(i+11,:)=data_R2_table(i,:);
		end
		Results_export(13+N,1)={'Kuiper V value'};
		Results_export(14+N,:)=colnames;
		for i=1:N
			Results_export(i+14+N,:)=data_V_table(i,:);
		end
		Results_export(16+2*N,1)={'KS D value'};
		Results_export(17+2*N,:)=colnames;
		for i=1:N
			Results_export(i+17+2*N,:)=data_D_table(i,:);
		end
		
		% Plot all results against mixed sample
		for i = 1:1
			cdf_D_best(:,i) = sum(cdf_source.*repmat(D_conc_out(i,2:end).*0.01,length(x_all),1),2);
		end
		axes(handles.KS_plot);
		hold on
		p2 = plot(x_all,cdf_sink, 'k', 'linewidth', 2);
		for i = 1:1
			p1 = plot(x_all,cdf_D_best(:,i), 'color', [0 1 0],'linewidth',1.5);
		end
		hold off
		axis([PDP_min PDP_max 0 1])
		title('KS test D statistic')
		xlabel('Age (Ma)')
		ylabel('Cumulative probability');
		legend([p1 p2],{'Optimized Model','Mixed sample'})
		
		for i = 1:1
			cdf_V_best(:,i) = sum(cdf_source.*repmat(V_conc_out(i,2:end).*0.01,length(x_all),1),2);
		end
		axes(handles.Kuiper_plot);
		hold on
		p4 = plot(x_all,cdf_sink, 'k', 'linewidth', 2);
		for i = 1:1
			p3 = plot(x_all,cdf_V_best(:,i),'color', [0 1 0],'linewidth',1.5);
		end
		hold off
		axis([PDP_min PDP_max 0 1])
		title('Kuiper test V statistic')
		xlabel('Age (Ma)')
		ylabel('Cumulative probability');
		legend([p3 p4],{'Optimized Model','Mixed sample'})
		
		for i = 1:1
			pdp_R2_best(:,i) = sum(pdp_source.*repmat(R2_conc_out(i,2:end).*0.01,length(x),1),2);
		end
		axes(handles.R2_plot);
		hold on
		p6 = plot(x,pdp_sink, 'k', 'linewidth', 2);
		for i = 1:1
			p5 = plot(x,pdp_R2_best(:,i),  'color', [0 1 0],'linewidth',1.5);
		end
		hold off
		axis([PDP_min PDP_max 0 max(pdp_sink)+max(pdp_sink)*.25])
		title('Cross-correlation coefficient')
		xlabel('Age (Ma)')
		ylabel('Relative probability');
		legend([p5 p6],{'Optimized Model','Mixed sample'})
		
		%update handles
		handles.D_conc_out=D_conc_out;
		handles.V_conc_out=V_conc_out;
		handles.R2_conc_out=R2_conc_out;
		handles.weights_OPT_best = best_OPTIMIZE1_wghts;
		
		% End of Optimization routine #1
		
		%% OPTION TO OPTIMIZE BASED ON BEST FIT MONTE CARLO METHOD UNMIXING RESULTS %%
	case handles.Matlab_optimize
		set(handles.Recent_run, 'String', 'Minimum function search optimization model');
		% Get best fits (mean values) from Monte Carlo results. These will be 'initial guesses' in the constrained minimum find algorithm below (fmincon).
		wghts_mc_D_tmp = sortrows([Passed_D, Passed_wghts_D],1);
		wghts_mc_D = wghts_mc_D_tmp(1:Max_best_fits,2:N+1).*100;
		wghts_mc_V_tmp = sortrows([Passed_V, Passed_wghts_V],1);
		wghts_mc_V = wghts_mc_V_tmp(1:Max_best_fits,2:N+1).*100;
		wghts_mc_R2_tmp = flipud(sortrows([Passed_R2, Passed_wghts_R2],1));
		wghts_mc_R2 = wghts_mc_R2_tmp(1:Max_best_fits,2:N+1).*100;
		
		lb = zeros(1,N); % Set lower bound (0%) for possible source contribution
		ub = ones(1,N)*100; % Set upper bound (100%) for possible source contribution
		opts1=  optimset('display','off');
		
		count = 0;
		h=waitbar(0,'Minimum Search Function','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
		setappdata(h,'canceling',0);
		% Apply constrained minimum algorithm. Each metric (D, V, and R2) calls a different function file (e.g., fmincon_D for KS D) that will use initial guesses (best mean values from Monte Carlo results).
		for i = 1:Max_best_fits
			count = count + 1;
			% Check for Cancel button press
			if getappdata(h,'canceling')
				break
			end
			[aa_D, bb_D] = fmincon(@fmincon_D, wghts_mc_D(i,:), [], [], [], [], lb, ub, @nonlcon2, opts1);
			ba_conc_D(i,:) = [bb_D,aa_D];
			
			[aa_V, bb_V] = fmincon(@fmincon_V, wghts_mc_V(i,:), [], [], [], [], lb, ub, @nonlcon2, opts1);
			ba_conc_V(i,:) = [bb_V,aa_V];
			
			[aa_R2, bb_R2] = fmincon(@fmincon_R2, wghts_mc_R2(i,:), [], [], [], [], lb, ub, @nonlcon2, opts1);
			ba_conc_R2(i,:) = [1-bb_R2,aa_R2];
			
			percent_complete = round((count)/Max_best_fits*100);
			waitbar(percent_complete/100,h);
		end
		delete(h)
		% Concatenate results
		ba_sort_D = sortrows(ba_conc_D,1);
		ba_sort_V = sortrows(ba_conc_V,1);
		ba_sort_R2 = flipud(sortrows(ba_conc_R2,1));
		assignin('base','ba_sort_D',ba_sort_D)
		
		% Sort and trim results
		D_mean_OPTIMIZE2 = mean(ba_sort_D(1:1,1));
		D_std_OPTIMIZE2 = std(ba_sort_D(1:1,1));
		V_mean_OPTIMIZE2 = mean(ba_sort_V(1:1,1));
		V_std_OPTIMIZE2 = std(ba_sort_V(1:1,1));
		R2_mean_OPTIMIZE2 = mean(ba_sort_R2(1:1,1));
		R2_std_OPTIMIZE2 = std(ba_sort_R2(1:1,1));
		assignin('base','D_mean_OPTIMIZE2',D_mean_OPTIMIZE2);
		
		% Calculate mean and stdev. of best fits and best overall fit. Print results in MATLAB command window
		mean_stdev_OPTIMIZE2_D_V_R2 = [[D_mean_OPTIMIZE2; D_std_OPTIMIZE2], [V_mean_OPTIMIZE2; V_std_OPTIMIZE2], [R2_mean_OPTIMIZE2; R2_std_OPTIMIZE2]];
		best_OPTIMIZE2_D_V_R2 = [ba_sort_D(1,1), ba_sort_V(1,1), ba_sort_R2(1,1)];
		best_OPTIMIZE2_wghts = [ba_sort_D(1,2:end); ba_sort_V(1,2:end); ba_sort_R2(1,2:end)].*.01;
		
		set(handles.mean_R2, 'String', round(R2_mean_OPTIMIZE2*1000)/1000);
		set(handles.mean_V,'String',round(V_mean_OPTIMIZE2*1000)/1000);
		set(handles.mean_D,'String',round(D_mean_OPTIMIZE2*1000)/1000);
		
		% Plot all results
		for i = 1:1
			cdf_D_best(:,i) = sum(cdf_source.*repmat(ba_sort_D(i,2:end).*0.01,length(x_all),1),2);
		end
		axes(handles.KS_plot);
		hold on
		p2 = plot(x_all,cdf_sink, 'k', 'linewidth', 2);
		for i = 1:1
			p1 = plot(x_all,cdf_D_best(:,i), 'color', [0 1 0],'linewidth',1.5);
		end
		hold off
		axis([PDP_min PDP_max 0 1])
		title('KS test D statistic')
		xlabel('Age (Ma)')
		ylabel('Cumulative probability');
		legend([p1 p2],{'Optimized Models','Mixed sample'})
		
		for i = 1:1
			cdf_V_best(:,i) = sum(cdf_source.*repmat(ba_sort_V(i,2:end).*0.01,length(x_all),1),2);
		end
		axes(handles.Kuiper_plot);
		hold on
		p4 = plot(x_all,cdf_sink, 'k', 'linewidth', 2);
		for i = 1:1
			p3 = plot(x_all,cdf_V_best(:,i), 'color', [0 1 0],'linewidth',1.5);
		end
		hold off
		axis([PDP_min PDP_max 0 1])
		title('Kuiper test V statistic')
		xlabel('Age (Ma)')
		ylabel('Cumulative probability');
		legend([p3 p4],{'Optimized Models','Mixed sample'})
		
		for i = 1:1
			pdp_R2_best(:,i) = sum(pdp_source.*repmat(ba_sort_R2(i,2:end).*0.01,length(x),1),2);
		end
		axes(handles.R2_plot);
		hold on
		p6 = plot(x,pdp_sink, 'k', 'linewidth', 3);
		for i = 1:1
			p5 = plot(x,pdp_R2_best(:,i),'color', [0 1 0],'linewidth',1.5);
		end
		hold off
		axis([PDP_min PDP_max 0 max(pdp_sink)+max(pdp_sink)*.25])
		title('Cross-correlation coefficient')
		xlabel('Age (Ma)')
		ylabel('Relative probability');
		legend([p5 p6],{'Optimized Models','Mixed sample'})
		
		%output average and standard deviations of model contributions to tables in
		%GUI interface
		colnames={'Samples', 'Relative Percent Contribution', 'Standard Deviation'};
		rows=transpose(headers);
		rows=rows(2:end,1);
		
		%if length(ba_sort_D) >= Max_best_fits
		Results_D = ba_sort_D(1,2:end);
		assignin('base','Results_D',Results_D)
		Results_D_mean_std = transpose([Results_D;zeros(1,N)]);
		assignin('base','Results_D_mean_std',Results_D_mean_std)
		data_D=num2cell(Results_D_mean_std);
		data_D_table=horzcat(rows,data_D);
		assignin('base','data_D_table',data_D_table)
		set(handles.D_table,'data',data_D_table, 'ColumnName', colnames)
		%end
		%if length(ba_sort_V) >= Max_best_fits
		Results_V = ba_sort_V(1,2:end);
		Results_V_mean_std = transpose([Results_V;zeros(1,N)]);
		data_V=num2cell(Results_V_mean_std);
		data_V_table=horzcat(rows,data_V);
		set(handles.V_table,'data',data_V_table, 'ColumnName', colnames)
		%end
		%if length(ba_sort_R2) >= Max_best_fits
		Results_R2 = ba_sort_R2(1,2:end);
		Results_R2_mean_std = transpose([Results_R2;zeros(1,N)]);
		data_R2=num2cell(Results_R2_mean_std);
		data_R2_table=horzcat(rows,data_R2);
		set(handles.R2_table,'data',data_R2_table, 'ColumnName', colnames)
		%end
		
		%Concatenate results for export
		Results_export= cell(14+3*N,3);
		Results_export(1,1)={'Results of Matlab Optimization unmixing model'};
		%if Licht_approach==0
		rad_on_source_scaling=get(handles.subsample_ages,'Value');
		if rad_on_source_scaling==0
			Results_export(2,1)={'weightings applied to PDPs and CDFs'};
			Results_export(3,:)={'number of grains in each trial=' 'N/A' ''};
		else
			
			age_num = str2num(get(handles.num_ages,'String'));
			Results_export(2,1)={'weightings applied to raw ages'};
			Results_export(3,:)={'number of grains in each trial=' age_num ''};
		end % End switch
		Results_export(4,:)={'trials=' trials ''};
		Results_export(5,:)={'percent accepted trials=' percentnum ''};
		Results_export(6,:)={'Cross-correlation cutoff=' threshold_R2 ''};
		Results_export(7,:)={'Kuiper V value cutoff=' threshold_V ''};
		Results_export(8,:)={'KS D value cutoff=' threshold_D ''};
		Results_export(10,:)={'Cross-correlation' '' ''};
		Results_export(11,:)=colnames;
		for i=1:N
			Results_export(i+11,:)=data_R2_table(i,:);
		end
		Results_export(13+N,1)={'Kuiper V value'};
		Results_export(14+N,:)=colnames;
		for i=1:N
			Results_export(i+14+N,:)=data_V_table(i,:);
		end
		Results_export(16+2*N,1)={'KS D value'};
		Results_export(17+2*N,:)=colnames;
		for i=1:N
			Results_export(i+17+2*N,:)=data_D_table(i,:);
		end
		handles.weights_OPT_best = best_OPTIMIZE2_wghts;
end % End if Optimization routine #2


set(handles.V_std,'String','');
set(handles.D_std,'String','');
set(handles.R2_std, 'String', '');
set(handles.D_std_dev2, 'String', 'Best D value:');
set(handles.V_std_dev2, 'String', 'Best V value:');
set(handles.R2_std_dev2, 'String', 'Best Cross-correlation coefficient:');
set(handles.plus_D, 'String', '');
set(handles.plus_V, 'String', '');
set(handles.plus_R2, 'String', '');

%update handles
%assignin('base','Results_R2',Results_R2)

handles.cdf_V_export=cdf_V_best;
handles.cdf_D_export=cdf_D_best;
handles.pdp_R2_export=pdp_R2_best;
handles.Results_D=Results_D;
handles.Results_V=Results_V;
handles.Results_R2=Results_R2;
handles.data_D=data_D;
handles.data_V=data_V;
handles.data_R2=data_R2;
handles.Results_export=Results_export;
handles.cdf_V_best=cdf_V_best;
handles.cdf_D_best=cdf_D_best;
handles.pdp_R2_best=pdp_R2_best;
handles.x_all=x_all;
guidata(hObject,handles);

% --- Executes on button press in Optimized_Licht.
function Optimized_Licht_Callback(hObject, eventdata, handles)
% hObject    handle to Optimized_Licht (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%cla(handles.R2_plot,'reset');
%cla(handles.KS_plot,'reset');
%cla(handles.Kuiper_plot,'reset');
set(handles.mean_V,'String','')
set(handles.mean_D,'String','')
set(handles.mean_R2,'String','')
set(handles.R2_std,'String','')
set(handles.V_std,'String','')
set(handles.D_std,'String','')




Licht_n=str2num(get(handles.num_ages,'String'));
Licht_N=str2num(get(handles.Licht_N,'String'));
trials = str2num(get(handles.run,'String'));
threshold_D = str2num(get(handles.V,'String'));
threshold_V = str2num(get(handles.V,'String'));
threshold_R2 = str2num(get(handles.R2_PDP,'String'));
PDP_min = str2num(get(handles.PDP_min,'String'));
PDP_max = str2num(get(handles.PDP_max,'String'));
PDP_step = str2num(get(handles.PDP_step,'String'));
x_axis_min = PDP_min;
x_axis_max = PDP_max;

N=handles.N;
weights_OPT_best = handles.weights_OPT_best;
data=handles.data;
x=handles.x;
pdp_sink=handles.pdp_sink;
cdf_V_best=handles.cdf_V_best;
cdf_D_best=handles.cdf_D_best;
pdp_R2_best=handles.pdp_R2_best;
x_all=handles.x_all;

%% Licht on optimization only

wait=waitbar(0,'Calculating','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wait,'canceling',0);

for r = 1:3
	
	%Check for Cancel button press
	if getappdata(wait,'canceling')
		break
	end
	
	% First determine how many ages correspond to each source based on rounding the randomly determined trial weight
	
	num_Ages = round(weights_OPT_best(r,:).*Licht_n);
	
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
		break
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
	
	% Combine mixed sample distribution and randomly generated model distribution to generate CDF
	%curves with equally binned x axes (required to calculate D and V)
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
	
	if r == 1
		count = 0;
		for p = 1:Licht_N
			count=count+1;
			%Check for Cancel button press
			if getappdata(wait,'canceling')
				break
			end
			D_tmp(1,p) = max([max(cdf_sink - cdf_source(:,p))],[max(cdf_source(:,p) - cdf_sink)]);
			waitbar(count/(3*Licht_N),wait);
		end
		D_mean(1,:) = mean(D_tmp(1,:));
		D_std(1,:) = std(D_tmp(1,:));
	end
	
	if r == 2
		count = 0;
		for p = 1:Licht_N
			count=count+1;
			%Check for Cancel button press
			if getappdata(wait,'canceling')
				break
			end
			V_tmp(1,p) = max(cdf_sink - cdf_source(:,p)) + max(cdf_source(:,p) - cdf_sink);
			waitbar((Licht_N/(3*Licht_N))+(count/(3*Licht_N)),wait);
		end
		V_mean(1,:) = mean(V_tmp(1,:));
		V_std(1,:) = std(V_tmp(1,:));
	end
	
	if r == 3
		for p = 1:Licht_N
			count=count+1;
			%Check for Cancel button press
			if getappdata(wait,'canceling')
				break
			end
			R2_tmp(1,p) = ((sum((pdp_sink - mean(pdp_sink)).*(pdp_source(:,p) - mean(pdp_source(:,p)))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp_source(:,p) - ...
				mean(pdp_source(:,p))).*(pdp_source(:,p) - mean(pdp_source(:,p))))))))*((sum((pdp_sink - mean(pdp_sink)).*(pdp_source(:,p) - mean(pdp_source(:,p)))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*...
				(pdp_sink - mean(pdp_sink))))*(sum((pdp_source(:,p) - mean(pdp_source(:,p))).*(pdp_source(:,p) - mean(pdp_source(:,p))))))));
			waitbar((2*Licht_N/(3*Licht_N))+(count/(3*Licht_N)),wait);
		end
		R2_mean(1,:) = mean(R2_tmp(1,:));
		R2_std(1,:) = std(R2_tmp(1,:));
	end
end % End once reached number of user specified trials
delete (wait)

R2_mean_out=R2_mean;
D_mean_out=D_mean;
V_mean_out=V_mean;
R2_std_out=D_std;
D_std_out=V_std;
V_std_out=R2_std;

set(handles.D_std_dev2, 'String', 'Mean D value:');
set(handles.V_std_dev2, 'String', 'Mean V value:');
set(handles.R2_std_dev2, 'String', 'Mean Cross-correlation coefficient:');
set(handles.plus_D, 'String', '+/-');
set(handles.plus_V, 'String', '+/-');
set(handles.plus_R2, 'String', '+/-');
set(handles.mean_R2, 'String', round(R2_mean_out*1000)/1000)
set(handles.mean_V,'String',round(V_mean_out*1000)/1000)
set(handles.mean_D,'String',round(D_mean_out*1000)/1000)
set(handles.R2_std, 'String', R2_std_out)
set(handles.V_std,'String',V_std_out)
set(handles.D_std,'String',D_std_out)

% Plot all results
axes(handles.KS_plot);
hold on
for i = 1:Licht_N
	p1 = plot(x_all,cdf_source(:,i), 'color', [1 .4 .2],'linewidth',1.5);
end
p2 = plot(x_all,cdf_sink, 'k', 'linewidth', 2);
hold off
axis([PDP_min PDP_max 0 1])
title('KS test D statistic')
xlabel('Age (Ma)')
ylabel('Cumulative probability');
legend([p1 p2],{'Subsampled Optimized Weighting','Mixed sample'},'Location','southeast')

axes(handles.Kuiper_plot);
hold on
for i = 1:Licht_N
	p3 = plot(x_all,cdf_source(:,i), 'color', [1 .4 .2],'linewidth',1.5);
end
p4 = plot(x_all,cdf_sink, 'k', 'linewidth', 2);
hold off
axis([PDP_min PDP_max 0 1])
title('Kuiper test V statistic')
xlabel('Age (Ma)')
ylabel('Cumulative probability');
legend([p3 p4],{'Subsampled Optimized Weighting','Mixed sample'},'Location','southeast')

axes(handles.R2_plot);
hold on
for i = 1:Licht_N
	p5 = plot(x,pdp_source(:,i),'color', [1 .4 .2],'linewidth',1.5);
end
p6 = plot(x,pdp_sink, 'k', 'linewidth', 3);
hold off
axis([PDP_min PDP_max 0 max(pdp_sink)+max(pdp_sink)*.25])
title('Cross-correlation coefficient')
xlabel('Age (Ma)')
ylabel('Relative probability');
legend([p5 p6],{'Subsampled Optimized Weighting','Mixed sample'})

set(handles.Recent_run,'String', 'Random ages selected from optimized weighting')
D_Licht_Results = sortrows([D_mean, D_std, weights_OPT_best(1,:)],1);
V_Licht_Results = sortrows([V_mean, V_std, weights_OPT_best(2,:)],1);
R2_Licht_Results = sortrows([R2_mean, R2_std, weights_OPT_best(3,:)],1);

guidata(hObject,handles);



% --- Executes on button press in clrplot.
function clrplot_Callback(hObject, eventdata, handles)
% hObject    handle to clrplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.R2_plot,'reset');
cla(handles.KS_plot,'reset');
cla(handles.Kuiper_plot,'reset');
set(handles.mean_V,'String','')
set(handles.mean_D,'String','')
set(handles.mean_R2,'String','')
set(handles.R2_std,'String','')
set(handles.V_std,'String','')
set(handles.D_std,'String','')
set(handles.Recent_run, 'String', 'Most recent run appears here');
set(handles.D_table,'data','')
set(handles.V_table,'data','')
set(handles.R2_table,'data','')
set(handles.D_std_dev2, 'String', '');
set(handles.V_std_dev2, 'String', '');
set(handles.R2_std_dev2, 'String', '');
set(handles.plus_D, 'String', '');
set(handles.plus_V, 'String', '');
set(handles.plus_R2, 'String', '');
set(handles.mean_R2, 'String', '')
set(handles.mean_V,'String','')
set(handles.mean_D,'String','')
set(handles.R2_std, 'String', '')
set(handles.V_std,'String','')
set(handles.D_std,'String','')


%%%%%%%%%%---------Plot Relative Contributions R2 Button---------%%%%%%%%%%
% --- Executes on button press in Plot_cont_R2.
function Plot_cont_R2_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_cont_R2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=cell2mat(handles.data_R2);
headers=handles.headers;
N=handles.N;
sample_number=1:N;
figure;
colours = colormap(jet((N)));
colorbar;
marker_size=10;
hold on
dx=0.3
for i=1:N
	e1=errorbar(sample_number(1,i),data(i,1),data(i,2),'o','MarkerSize',...
		marker_size,'MarkerEdgeColor','black','MarkerFaceColor',colours((i),:));
end
for i=1:N
	e2=text(sample_number(1,i)+dx,data(i,1),headers(1,1+i));
	set(e2,'Rotation',90)
end

title('Relative contributions from Cross-correlation coefficient');
xlabel('Sample number');
ylabel('Relative contribution');
set(gca,'XTick',[1:N]);
hold off

%%%%%%%%%%---------Plot Relative Contributions D Button---------%%%%%%%%%%
% --- Executes on button press in Plot_cont_D.
function Plot_cont_D_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_cont_D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=cell2mat(handles.data_D);
headers=handles.headers;
N=handles.N;
sample_number=1:N;
figure;
colours = colormap(jet((N)));
colorbar;
marker_size=10;
hold on
dx=0.3
for i=1:N
	e1=errorbar(sample_number(1,i),data(i,1),data(i,2),'o','MarkerSize',...
		marker_size,'MarkerEdgeColor','black','MarkerFaceColor',colours((i),:));
end
for i=1:N
	e2=text(sample_number(1,i)+dx,data(i,1),headers(1,1+i));
	set(e2,'Rotation',90)
end
title('Relative contributions from KS test D value');
xlabel('Sample number');
ylabel('Relative contribution');
set(gca,'XTick',[1:N]);
hold off

%%%%%%%%%%---------Plot Relative Contributions V Button---------%%%%%%%%%%
% --- Executes on button press in Plot_cont_V.
function Plot_cont_V_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_cont_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=cell2mat(handles.data_V);
headers=handles.headers;
N=handles.N;
sample_number=1:N;
figure;
colours = colormap(jet((N)));
colorbar;
marker_size=10;
hold on
dx=0.3;
for i=1:N
	e1=errorbar(sample_number(1,i),data(i,1),data(i,2),'o','MarkerSize',...
		marker_size,'MarkerEdgeColor','black','MarkerFaceColor',colours((i),:));
end
for i=1:N
	e2=text(sample_number(1,i)+dx,data(i,1),headers(1,1+i));
	set(e2,'Rotation',90)
end
title('Relative contributions from Kuiper test V value');
xlabel('Sample number');
ylabel('Relative contribution');
set(gca,'XTick',[1:N]);
hold off

% --- Executes on button press in Export_figs.
function Export_figs_Callback(hObject, eventdata, handles)
% hObject    handle to Export_figs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pdps=handles.pdps;
data=handles.data;
nsamples=handles.nsamples;
x=handles.x;
PDP_min = str2num(get(handles.PDP_min,'String'));
PDP_max = str2num(get(handles.PDP_max,'String'));
Results_D=handles.Results_D
Passed_cdf_D=handles.Passed_cdf_D
plotnum=handles.plotnum
N=handles.N

f = figure;
R2_plot=findall(handles.R2_plot,'type','axes');
ax.Units='normalized';
ax.Position = [0 0 1 1];
copyobj(R2_plot,f);
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
legend({'Models','Mixed sample'})

g=figure;
Kuiper_plot=handles.Kuiper_plot
ax.Units='normalized';
ax.Position = [0 0 1 1];
copyobj(Kuiper_plot,g);
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
legend({'Models','Mixed sample'})

h=figure;
KS_plot=handles.KS_plot
ax.Units='normalized';
ax.Position = [0 0 1 1];
copyobj(KS_plot,h);
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
legend({'Models','Mixed sample'})

i=figure
colours = colormap(jet((nsamples-1)));
for i = 2:nsamples;
	datai = data(:,i*2-1);
	datai =datai(isfinite(datai(:,1)),:);
	cdf(i) = cdfplot(datai);
	set(cdf(i),'color',colours((i-1),:),'linewidth',1.5);
	hold on;
	grid on;
	xlabel('');
	ylabel('');
end
q2=cdfplot(data(:,1))
set(q2,'color','k','linewidth',2)
axis([PDP_min, PDP_max, 0, 1]);
legend([q2],{'Mixed sample'})
title('Source and Mixed CDFs');
xlabel('Age (Ma)')
ylabel('Cumulative probability');
colorbar;
hold off

j=figure
hold on
colours = colormap(jet((nsamples-1)));
pdp_source = pdps(:,2:end);
pdp_sink = pdps(:,1);
for i = 1:nsamples-1;
	plot(x,pdp_source(:,i),'color',colours((i),:),'linewidth',1.5);
	grid on
end
q4=plot(x,pdp_sink, 'k', 'linewidth', 2)
legend([q4],{'Mixed sample'})
title('Source and Mixed PDPs')
xlabel('Age (Ma)')
ylabel('Relative probability');
ax.YAxis.Exponent = -2;
colorbar;
hold off


% --- Executes on button press in Export_Data.
function Export_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Export_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Results_export=handles.Results_export;
[file,path] = uiputfile('*.xls','Save file');
%xlswrite([path file], Results_export);
writetable(table(Results_export),[path file])

% --- Executes on button press in Export_PDP_CDF.
function Export_PDP_CDF_Callback(hObject, eventdata, handles)
% hObject    handle to Export_PDP_CDF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
N=handles.N
x_all=handles.x_all;
plotnum=handles.plotnum;
%cdf_D_export=handles.cdf_D_export;
%cdf_V_export=handles.cdf_V_export;
%pdp_R2_export=handles.pdp_R2_export;
Results_D=transpose(handles.Results_D);
Results_V=transpose(handles.Results_V);
Results_R2=transpose(handles.Results_R2);
PDP_min = str2num(get(handles.PDP_min,'String'));
PDP_max = str2num(get(handles.PDP_max,'String'));
PDP_step = str2num(get(handles.PDP_step,'String'));
x_pdp=transpose(PDP_min:PDP_step:PDP_max);
N = str2num(get(handles.numsamples,'string'));
Results_R2_out = [x_pdp,Results_R2(3+N:end,:)];
[file,path] = uiputfile('model_PDPs_R2.xls','Save file');
%xlswrite([path file], Results_R2);
writetable(table(Results_R2_out),[path file]);
Results_V_out = [x_all,Results_V(3+N:length(x_all)+2+N,:)];
[file,path] = uiputfile('model_CDFs_Kuiper.xls','Save file');
%xlswrite([path file], Results_V);
writetable(table(Results_V_out),[path file]);
Results_D_out = [x_all,Results_D(3+N:length(x_all)+2+N,:)];
[file,path] = uiputfile('model_CDFs_KS.xls','Save file');
%xlswrite([path file], Results_D);
writetable(table(Results_D_out),[path file]);

%%%%%%%%%%%%%%%%---------Example Data Set Button-----------%%%%%%%%%%%%%%%%
% --- Executes on button press in Example.
function Example_Callback(hObject, eventdata, handles)
% hObject    handle to Example (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Example_Data_Set;


function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of run as text
%        str2double(get(hObject,'String')) returns contents of run as a double


% --- Executes during object creation, after setting all properties.
function run_CreateFcn(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end



function percentnum_Callback(hObject, eventdata, handles)
% hObject    handle to percentnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of percentnum as text
%        str2double(get(hObject,'String')) returns contents of percentnum as a double


% --- Executes during object creation, after setting all properties.
function percentnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to percentnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in skittle_plots.
function skittle_plots_Callback(hObject, eventdata, handles)
% hObject    handle to skittle_plots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=handles.data;
nsamples=handles.nsamples;
text1=handles.text1;
N=handles.N;
Results_D=handles.Results_D;
optimization=handles.optimization;
Passed_D=handles.Passed_D;
Passed_wghts_D=handles.Passed_wghts_D;
Passed_V=handles.Passed_V;
Passed_wghts_V=handles.Passed_wghts_V;
Passed_R2=handles.Passed_R2;
Passed_wghts_R2=handles.Passed_wghts_R2;

%assignin('base' ,'Passed_D', Passed_D)
%assignin('base', 'N', N)
%assignin('base','Results_D', Results_D)
if optimization>0
	D_conc_out=handles.D_conc_out;
	R2_conc_out=handles.R2_conc_out;
	V_conc_out=handles.V_conc_out;
end

for i=1:nsamples
	headers(1,i)=text1(1,2*i-1);
end

%calculate 1 or 2 sigma
rad_on=get(handles.sigma,'selectedobject');
switch rad_on
	case handles.sigma1
		%
	case handles.sigma2
		for i=1:nsamples
			data(:,2*i)=data(:,2*i)./2;
		end
	otherwise
		set(handles.edit_radioselect,'string','');
end

%Get user inputs
trials = str2num(get(handles.run,'String'));
plotnum = trials;
threshold_D = 1;
threshold_V = 1;
threshold_R2 = 0;
PDP_min = str2num(get(handles.PDP_min,'String'));
PDP_max = str2num(get(handles.PDP_max,'String'));
PDP_step = str2num(get(handles.PDP_step,'String'));
x_axis_min = PDP_min;
x_axis_max = PDP_max;

%Plot visualiztion plots
count=0;
w=waitbar(0,'Plotting figure 1 of 3: Visualization of Monte Carlo mixture modeling using KS test D value',...
	'CreateCancelBtn',...
	'setappdata(gcbf,''canceling'',1)');
movegui(w,[600,200])
setappdata(w,'canceling',0);
KS_D=figure;
colours = colormap(jet((nsamples-1)));
colorbar;
hold on
chk_bx=get(handles.watch,'Value');
if chk_bx==0
	set(KS_D,'Visible','off')
end
if optimization==0
	while count<trials;
		% Check for Cancel button press
		if getappdata(w,'canceling')
			break
		end
		count=count+1;
		for i=1:N
			line(Passed_wghts_D(count,i),Passed_D(count,1),'LineStyle','none',...
				'Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
				'MarkerSize',6);
		end
		waitbar(single(count/trials),w)
	end
else
	for i = 1:N
		for j = 1:length(D_conc_out(:,1))
			% Check for Cancel button press
			if getappdata(w,'canceling')
				break
			end
			line(D_conc_out(j,i+1),D_conc_out(j,1),'LineStyle','none',...
				'Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
				'MarkerSize',6);
			count=count+1;
			waitbar(count/(length(D_conc_out(:,1)*N)),w);
		end
	end
	
	
end
figure(KS_D)
set(0, 'CurrentFigure', KS_D)
set(gca, 'Ydir', 'reverse');
set(KS_D,'Visible','on')
xlabel('Relative contribution');
ylabel('KS test D value');
if optimization==0
	title('Monte Carlo mixture modeling using KS test D value');
else
	title('Optimized mixture modeling using KS test D value');
end
hold off
set(0,'CurrentFigure',w)
delete (w)

count=0;
w=waitbar(0,'Plotting figure 2 of 3: Visualization of Monte Carlo mixture modeling using Kuiper test V value',...
	'CreateCancelBtn',...
	'setappdata(gcbf,''canceling'',1)');
movegui(w,[600,200])
setappdata(w,'canceling',0);
marker_size=25;
Kuiper_V=figure;
colours = colormap(jet((nsamples-1)));
colorbar;
hold on
if chk_bx==0
	set(Kuiper_V,'Visible','off')
end
if optimization==0
	while count<trials
		% Check for Cancel button press
		if getappdata(w,'canceling')
			break
		end
		count=count+1;
		
		for i=1:N
			line(Passed_wghts_V(count,i),Passed_V(count,1),'LineStyle','none',...
				'Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
				'MarkerSize',6);
		end
		waitbar(single(count/trials),w)
	end
else
	for i = 1:N
		for j = 1:length(V_conc_out(:,1))
			% Check for Cancel button press
			if getappdata(w,'canceling')
				break
			end
			line(V_conc_out(j,i+1),V_conc_out(j,1),'LineStyle','none',...
				'Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
				'MarkerSize',6);
			count=count+1;
			waitbar(count/(length(D_conc_out(:,1)*N)),w);
		end
	end
end
figure(Kuiper_V)
set(0, 'CurrentFigure', Kuiper_V)
set(gca, 'Ydir', 'reverse');
set(Kuiper_V,'Visible','on')
xlabel('Relative contribution');
ylabel('Kuiper test V value');
if optimization==0
	title('Monte Carlo mixture modeling using Kuiper test V value');
else
	title('Optimized mixture modeling using Kuiper test V value');
end
hold off
set(0,'CurrentFigure',w)
delete (w)

count=0;
w=waitbar(0,'Plotting figure 3 of 3: Vizualization of Monte Carlo mixture modeling using Cross-correlation coefficient',...
	'CreateCancelBtn',...
	'setappdata(gcbf,''canceling'',1)');
movegui(w,[600,200])
setappdata(w,'canceling',0);
marker_size=25;
R2=figure;
colours = colormap(jet((nsamples-1)));
colorbar;
hold on
if chk_bx==0
	set(R2,'Visible','off')
end
if optimization==0
	while count<trials
		% Check for Cancel button press
		if getappdata(w,'canceling')
			break
		end
		count=count+1;
		
		for i=1:N
			line(Passed_wghts_R2(count,i),Passed_R2(count,1),'LineStyle','none',...
				'Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
				'MarkerSize',6);
		end
		waitbar(single(count/trials),w)
	end
else
	for i = 1:N
		for j = 1:length(R2_conc_out(:,1))
			% Check for Cancel button press
			if getappdata(w,'canceling')
				break
			end
			line(R2_conc_out(j,i+1),R2_conc_out(j,1),'LineStyle','none',...
				'Marker','o','MarkerFaceColor',colours((i),:),'MarkerEdgeColor','none',...
				'MarkerSize',6);
			count=count+1;
			waitbar(count/(length(D_conc_out(:,1)*N)),w);
		end
	end
	
end
figure(R2)
set(R2,'Visible','on')
xlabel('Relative contribution');
ylabel('Cross-correlation coefficient');
if optimization==0
	title('Monte Carlo mixture modeling using Cross-correlation coefficient');
else
	title('Optimized mixture modeling using Cross-correlation coefficient');
end
hold off
set(0,'CurrentFigure',w)
delete (w)


% --- Executes on button press in Raw_ages.
function Raw_ages_Callback(hObject, eventdata, handles)
% hObject    handle to Raw_ages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Raw_ages



function num_ages_Callback(hObject, eventdata, handles)
% hObject    handle to num_ages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_ages as text
%        str2double(get(hObject,'String')) returns contents of num_ages as a double


% --- Executes during object creation, after setting all properties.
function num_ages_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_ages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end



function R2_PDP_Callback(hObject, eventdata, handles)
% hObject    handle to R2_PDP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R2_PDP as text
%        str2double(get(hObject,'String')) returns contents of R2_PDP as a double


% --- Executes during object creation, after setting all properties.
function R2_PDP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R2_PDP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end



function V_Callback(hObject, eventdata, handles)
% hObject    handle to V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of V as text
%        str2double(get(hObject,'String')) returns contents of V as a double


% --- Executes during object creation, after setting all properties.
function V_CreateFcn(hObject, eventdata, handles)
% hObject    handle to V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end



function D_Callback(hObject, eventdata, handles)
% hObject    handle to D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of D as text
%        str2double(get(hObject,'String')) returns contents of D as a double


% --- Executes during object creation, after setting all properties.
function D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end



function PDP_step_Callback(hObject, eventdata, handles)
% hObject    handle to PDP_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PDP_step as text
%        str2double(get(hObject,'String')) returns contents of PDP_step as a double


% --- Executes during object creation, after setting all properties.
function PDP_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PDP_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end



function PDP_max_Callback(hObject, eventdata, handles)
% hObject    handle to PDP_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PDP_max as text
%        str2double(get(hObject,'String')) returns contents of PDP_max as a double


% --- Executes during object creation, after setting all properties.
function PDP_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PDP_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end



function PDP_min_Callback(hObject, eventdata, handles)
% hObject    handle to PDP_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PDP_min as text
%        str2double(get(hObject,'String')) returns contents of PDP_min as a double


% --- Executes during object creation, after setting all properties.
function PDP_min_CreateFcn(hObject, ~, handles)
% hObject    handle to PDP_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in watch.
function watch_Callback(hObject, eventdata, handles)
% hObject    handle to watch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of watch

function Max_best_fits_Callback(hObject, eventdata, handles)
% hObject    handle to Max_best_fits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Max_best_fits as text
%        str2double(get(hObject,'String')) returns contents of Max_best_fits as a double


% --- Executes during object creation, after setting all properties.
function Max_best_fits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Max_best_fits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

function smallest_text_Callback(hObject, eventdata, handles)
% hObject    handle to Recent_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smallest_text as text
%        str2double(get(hObject,'String')) returns contents of smallest_text as a double

% --- Executes during object creation, after setting all properties.
function smallest_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smallest_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

function Recent_run_Callback(hObject, eventdata, handles)
% hObject    handle to Recent_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Recent_run as text
%        str2double(get(hObject,'String')) returns contents of Recent_run as a double


% --- Executes during object creation, after setting all properties.
function Recent_run_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Recent_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end



function Licht_N_Callback(hObject, eventdata, handles)
% hObject    handle to Licht_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Licht_N as text
%        str2double(get(hObject,'String')) returns contents of Licht_N as a double


% --- Executes during object creation, after setting all properties.
function Licht_N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Licht_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Optimized_Licht.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to Optimized_Licht (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in skittle_plots.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to skittle_plots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in subsample_ages.
function subsample_ages_Callback(hObject, eventdata, handles)
% hObject    handle to subsample_ages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subsample_ages



function kde_bandwidth_Callback(hObject, eventdata, handles)
% hObject    handle to kde_bandwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kde_bandwidth as text
%        str2double(get(hObject,'String')) returns contents of kde_bandwidth as a double


% --- Executes during object creation, after setting all properties.
function kde_bandwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kde_bandwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end
