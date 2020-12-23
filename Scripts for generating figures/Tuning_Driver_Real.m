% LayerOptimizer Tuning
%
% Jon Cooper

set(0, 'DefaultLineLineWidth', 2);
% close all
rng(0);

%% Load dataset

filename = 'cpt-data-Netherlands/CPT_45185_Raw01.TXT';
FID=fopen(char(filename));
datacell = textscan(FID, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%q%q%f%q%q%f%f%f%q%f',...
    'HeaderLines', 17);
fclose(FID);
%
qcToAtm = 1000.0/101.3;
fsToAtm = 1000.0/101.3;
Patm    = 101.3;
gamma1  = 18.0;
gamma2  = 19.0;
gammaw  = 9.81;
zwater  = 1.0;
z       = datacell{1};
qcMeas  = datacell{9}*qcToAtm; % THIS ONE
fs      = datacell{3}*fsToAtm;
%
% ===================================================================
%  Impose limits on qc and fs to avoid any negative values
%       - including at the end of the fs file
%
qcMeas  = max(0.01,qcMeas);
fs      = max(0.001,fs);
qlength = length(qcMeas);
dz      = (max(z)-min(z))/(length(z)-1);
%
% ===================================================================
%  Inputs for analyses & filtering
%
coneD  = 0.03568;
dzd    = dz/coneD;
zd     = z/coneD;
obs = [zd,qcMeas];

trunct = 50; truncb = 0;
obs = obs(trunct:end-truncb,:); zd = zd(trunct:end-truncb,:);% qcTrue = qcTrue(trunct:end-truncb); %qcMeas = qcMeas(50:end-50);
zd = zd-min(zd);

MAX = 450;

% Clear temporary variables
clear opts SoilModel IA zdt zdtt zdm

%% Artificial True Model (sine wave, 3 different resistances)
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
% qcMeas = blur(qcTrue);
% obs(:,2) = qcMeas;

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

options = {'MaxIter',10;'MIN',.01;'MAX',MAX;'AOIMaxAdded',2;'LayerTol',0.003;'Plt',1};
[layers,info] = LayerOptimizerLog(obs,blur,options);

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
plot(obs(:,2),-zd,'--b'), %plot(qcTrue,-zd,'-b'), hold off;
legend('Inverse Model','Top of Layer','Blur of Model','Observed','True','location','southwest');
title("Original Model");
ylabel('Depth (m)','interpreter','LaTeX'); xlabel('$q_c$ Resistance (MPa)','interpreter','LaTeX');
xlim([0,MAX]); ylim([-max(zd),0]);
set(gca,'FontSize',12);
subplot(1,2,2);
t = LayerModelEval(layers,zd);
plot(t,-zd,'-r'), hold on, plot(layers(Nf+1:2*Nf),-layers(1:Nf),'or'), plot(blur(t),-zd,'--r');
plot(obs(:,2),-zd,'--b'), %plot(qcTrue,-zd,'-b'), hold off;
legend('Inverse Model','Top of Layer','Blur of Model','Observed','True','location','southwest');
title("Log Model");
ylabel('Depth (m)','interpreter','LaTeX'); xlabel('$q_c$ Resistance (MPa)','interpreter','LaTeX');
xlim([0,MAX]); ylim([-max(zd),0]);
set(gca,'FontSize',12);

% subplot(1,3,3);
% plot(t,-zd,'-r'), hold on, plot(layers(Nf+1:2*Nf),-layers(1:Nf),'or'), plot(mask*conv(t,kern),-zd,'--r');
% plot(obs(:,2),-zd,'--b'), plot(qcTrue,-zd,'-b'), hold off;
% legend('Inverse Model','Top of Layer','Blur of Model','Observed','True','location','southwest');
% title(Nf+" Layers, Final Model");
% ylabel('Depth (m)','interpreter','LaTeX'); xlabel('$q_c$ Resistance (MPa)','interpreter','LaTeX');
% xlim([0,MAX]);
% set(gca,'FontSize',12);

