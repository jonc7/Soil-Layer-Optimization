% LayerOptimizer Tuning
%
% Jon Cooper

set(0, 'DefaultLineLineWidth', 2);
% close all
rng(0);

%% Load dataset
MODEL = 4; CPT = 2; %9;3 for hard
% DATA = 7; SHEET = 1;
% load('CPT40True.mat');
% obs = [CPT40True.zd,CPT40True.qcMeas]; % observed soil resistances
% trueProfile = CPT40True.qcTrue; % true resistance profile
MAX = 60; MIN = .01;

% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 2);
opts.Sheet = "SM "+MODEL+"- CPT"+CPT;
opts.DataRange = "A4:B700";
opts.VariableNames = ["A", "B"];
opts.SelectedVariableNames = ["A", "B"];
opts.VariableTypes = ["double", "double"];

% Import the data
SoilModel = readtable("C:\Users\retin\Documents\Thin Layer Inversion\cpt-data-Netherlands\TrueData_GoodMatches.xlsx", opts, "UseExcel", false);
SoilModel = SoilModel{:,:};
zdm = SoilModel(:,1); qcMeas = SoilModel(:,2);
zdm = zdm(~isnan(qcMeas)); qcMeas = qcMeas(~isnan(qcMeas));
zdm = zdm(~isnan(zdm)); qcMeas = qcMeas(~isnan(zdm));

% Repeat for True data
opts = spreadsheetImportOptions("NumVariables", 2);
opts.Sheet = "SM "+MODEL+"- CPT"+CPT;
opts.DataRange = "G4:H700";
opts.VariableNames = ["G", "H"];
opts.SelectedVariableNames = ["G", "H"];
opts.VariableTypes = ["double", "double"];
SoilModel = readtable("C:\Users\retin\Documents\Thin Layer Inversion\cpt-data-Netherlands\TrueData_GoodMatches.xlsx", opts, "UseExcel", false);
SoilModel = SoilModel{:,:};
zdt = SoilModel(:,1); qcTrue = SoilModel(:,2);
zdt = zdt(~isnan(qcTrue)); qcTrue = qcTrue(~isnan(qcTrue));
zdt = zdt(~isnan(zdt)); qcTrue = qcTrue(~isnan(zdt));
[zdtt,IA,~] = unique(zdt); zdt = zdt+0.0001; zdt(IA) = zdtt-0.0001; % handle duplicate layers
[zdt,IA,~] = unique(zdt); qcTrue = qcTrue(IA);

% Interpolate to match length of measured data & sampling frequency
if length(zdm) > length(zdt)
    M = length(zdm);
    zd = linspace(min(zdm),max(zdm),M)';
else
    M = length(zdt);
    zd = linspace(min(zdt),max(zdt),M)';
end

qcTrue = interp1(zdt,qcTrue,zd,'linear','extrap'); %zd = zd(~isnan(qcMeas)); qcTrue = qcTrue(~isnan(qcTrue));
qcMeas = interp1(zdm,qcMeas,zd,'linear','extrap'); %qcMeas = qcMeas(~isnan(qcMeas));
% qcTrue = qcTrue/100; qcMeas = qcMeas/100;
obs = [zd,qcMeas];

trunct = 50; truncb = 0;
obs = obs(trunct:end-truncb,:); zd = zd(trunct:end-truncb,:); qcTrue = qcTrue(trunct:end-truncb); %qcMeas = qcMeas(50:end-50);
zd = zd-min(zd);

% Clear temporary variables
clear opts SoilModel IA zdt zdtt zdm

%% Random True Models
% % N = floor(10+11*rand(1));
% N = 15;
% % pos = max(zd)*sort(rand(N,1));
% pos = max(zd)*(0:N-1)'/N/2;
% % res = 5+55*rand(N,1);
% res = [10;60]; res = [res;res;res;res;res;res;res;res;res;res]; res = res(1:N);
% qcTrue = LayerModelEval([pos;res],zd);
% plot(qcTrue,-zd);

%% Artificial True Model
% Ti = (1/30/pi); Tf = (pi/4);
% pos = Ti:((Tf-Ti)/(length(zd)-1)):Tf;
% res = feval(@(x) sin(1./x),pos);
% pos = pos(nonzeros((1:length(pos)-1).*abs((res(1:end-1) >= 0) - (res(2:end) >= 0))+...
%     +(1:length(pos)-1).*abs((res(1:end-1) < 0) - (res(2:end) < 0))))';
% res = repmat(MAX*[.1;.45],ceil(length(pos)/2),1);
% res = res(1:length(pos));
% 
% pos2 = Ti:((Tf-Ti)/(length(zd)-1)):Tf;
% res2 = feval(@(x) cos(1./x),pos2);
% pos2 = pos2(nonzeros((1:length(pos2)-1).*abs((res2(1:end-1) > 0) - (res2(2:end) > 0))+...
%     +(1:length(pos2)-1).*abs((res2(1:end-1) < 0) - (res2(2:end) < 0)))/2)';
% res2 = repmat(MAX*[.1;.9],ceil(length(pos2)/2),1);
% res2 = res2(1:length(pos2));
% 
% qcTrue = LayerModelEval([pos;pos2;res;res2],zd);

%% PSF Setup
% kern = normpdf(linspace(-3,3,65)); % use a basic bell curve for the psf
kern = chi2pdf(linspace(0,8,65),4);
kern = kern/sum(kern);

%% Force Blur Measures
M = size(obs,1); mask = eye(M); n = 0; z = 0;
mask = [zeros(M,floor(length(kern)/2)),mask,zeros(M,ceil(length(kern)/2)-1)];
mask = diag([zeros(1,z),linspace(0,1,n),ones(1,M-n-z)].*[ones(1,M-n-z),linspace(1,0,n),zeros(1,z)])*mask;
blur = @(layer) mask*conv(layer,kern);
qcMeas = blur(qcTrue);
obs(:,2) = qcMeas;

%% Layer0 Automatic Selection
% pos = (qcMeas(3:end)-qcMeas(2:end-1) >= 0) & (qcMeas(1:end-2)-qcMeas(2:end-1) > 0);
% neg = (qcMeas(3:end)-qcMeas(2:end-1) <= 0) & (qcMeas(1:end-2)-qcMeas(2:end-1) < 0);
% loc = nonzeros((pos | neg).*zd(2:end-1)); res = nonzeros((pos| neg).*qcMeas(2:end-1));
% loc = loc - (zd(length(kern))-zd(1))/4;
% loc = [min(zd);loc]; res = [qcMeas(ceil(length(kern)/2));res];
% N = length(loc);

% manually input layers
% loc(end+1) = .07; loc(end+1) = .5;
% res(end+1) = 80; res(end+1) = 80;

% layers0 = [pos;res];

% NOTE: be sure to include ';'layers0',layers0;' in the options!
% options = [options;{'layers0',layers}];

%% Optimize

options = {'MaxIter',4;'MIN',.04;'MAX',MAX;'Plt',0;'parallel',0};%;'layers0',layers0};
tic; [layers,info] = LayerOptimizerLog(obs,blur,options);

%% Plot
layers0 = info.layers0;
N0 = info.Ninitial;
Nf = info.Nfinal;
IM = info.InitialMisfit;
FM = info.FinalMisfit;

fprintf('Initial Misfit: %f\n',IM);
fprintf('Final Misfit: %f\n',FM);

figure
subplot(1,2,1);
t0 = LayerModelEval(layers0,zd);
plot(t0,-zd,'-r'), hold on, plot(layers0(N0+1:2*N0),-layers0(1:N0),'or'), plot(blur(t0),-zd,'--r');
plot(obs(:,2),-zd,'--b'), plot(qcTrue,-zd,'-b'), hold off;
legend('q_c^{inv}','Top of Layer','d^{sim}(q_c^{inv})','d^{meas}','q_c^{true}','location','southwest');
title("Initial Guess");
ylabel('Depth (m)','interpreter','LaTeX'); xlabel('$q_{c1n}$ Resistance','interpreter','LaTeX');
xlim([0,MAX]); ylim([-max(zd),0]);
set(gca,'FontSize',12);
subplot(1,2,2);
t = LayerModelEval(layers,zd);
plot(t,-zd,'-r'), hold on, plot(layers(Nf+1:2*Nf),-layers(1:Nf),'or'), plot(blur(t),-zd,'--r');
plot(obs(:,2),-zd,'--b'), plot(qcTrue,-zd,'-b'), hold off;
legend('q_c^{inv}','Top of Layer','d^{sim}(q_c^{inv})','d^{meas}','q_c^{true}','location','southwest');
title("Log Misfit");
ylabel('Depth (m)','interpreter','LaTeX'); xlabel('$q_{c1n}$ Resistance','interpreter','LaTeX');
xlim([0,MAX]); ylim([-max(zd),0]);
set(gca,'FontSize',12);

