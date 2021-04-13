% LayerOptimizerLog_Driver
% example usage of LayerOptimizerLog.m
%
% Jon Cooper

set(0, 'DefaultLineLineWidth', 2);
close all
rng(0);

%% Load dataset
load('cpt-data-Netherlands\CPT40True.mat'); % to be replaced
obs = [CPT40True.zd,CPT40True.qcMeas]; % observed soil resistances
trueProfile = CPT40True.qcTrue; % true resistance profile
M = length(trueProfile);

%% Setup
% kern = normpdf(linspace(-3,3,65)); % use a basic bell curve for the psf
kern = chi2pdf(linspace(0,8,65),4); % use chi^2 curve for asymmetry
kern = kern/sum(kern);
MAX = 5;
options = {'ConvTol',1e-3;'Plt',1;'MAX',MAX};

% !temporary! force blur to match psf type for demo purposes only
mask = eye(size(obs,1)); n = 0; z = 0;
mask = [zeros(size(obs,1),floor(length(kern)/2)),mask,zeros(size(obs,1),ceil(length(kern)/2)-1)];
mask = diag([zeros(1,z),linspace(0,1,n),ones(1,M-n-z)].*[ones(1,M-n-z),linspace(1,0,n),zeros(1,z)])*mask;
blur = @(layers) mask*conv(layers,kern);
obs(:,2) = blur(trueProfile);

%% Optimize
[layers,info] = LayerOptimizerLog(obs,blur,options);

%% Plot
zd = CPT40True.zd;

layers0 = info.layers0;
Ni = info.Ninitial;
Nf = info.Nfinal;
IM = info.InitialMisfit;
FM = info.FinalMisfit;

fprintf('Initial Misfit: %f\n',IM);
fprintf('Final Misfit: %f\n',FM);

figure
subplot(1,2,1);
t0 = LayerModelEval(layers0,zd);
plot(t0,-zd,'-r'), hold on, plot(layers0(Ni+1:2*Ni),-layers0(1:Ni),'or'), plot(blur(t0),-zd,'--r');
plot(obs(:,2),-zd,'--b'), plot(trueProfile,-zd,'-b'), hold off;
legend('$q_c$','Top of Layer','$\tilde{q_c}^{sim}$','$\tilde{q_c}^{meas}$','$q_c^{true}$','interpreter','latex');
title('Initial Guess','interpreter','latex');
xlim([0,MAX]);
subplot(1,2,2);
t = LayerModelEval(layers,zd);
plot(t,-zd,'-r'), hold on, plot(layers(Nf+1:2*Nf),-layers(1:Nf),'or'), plot(blur(t),-zd,'--r');
plot(obs(:,2),-zd,'--b'), plot(trueProfile,-zd,'-b'), hold off;
legend('$q_c$','Top of Layer','$\tilde{q_c}^{sim}$','$\tilde{q_c}^{meas}$','$q_c^{true}$','interpreter','latex');
title('Final $q_c$','interpreter','latex');
xlim([0,MAX]);
