% LayerOptimizer_Driver
% example usage of LayerOptimizer.m
%
% Jon Cooper

set(0, 'DefaultLineLineWidth', 2);
close all
rng(0);

%% Load dataset
load('CPT40True.mat');
obs = [CPT40True.zd,CPT40True.qcMeas]; % observed soil resistances
trueProfile = CPT40True.qcTrue; % true resistance profile

%% Setup
kern = normpdf(linspace(-3,3,65)); % use a basic bell curve for the psf
kern = kern/sum(kern);
N = 14; % initial number of layers
options = {'ConvTol',1e-3;'Plt',0};

% !temporary! force blur to match psf type for demo purposes only
mask = eye(size(obs,1));
mask = [zeros(size(obs,1),floor(length(kern)/2)),mask,zeros(size(obs,1),ceil(length(kern)/2)-1)];
obs(:,2) = mask*conv(trueProfile,kern);

%% Optimize
[layers,info] = LayerOptimizer(N,obs,kern,options);

%% Plot
zd = CPT40True.zd;

layers0 = info.layers0;
Nf = info.Nfinal;
IM = info.InitialMisfit;
FM = info.FinalMisfit;
mask = info.mask;

fprintf('Initial Misfit: %f\n',IM);
fprintf('Final Misfit: %f\n',FM);

figure
subplot(1,2,1);
t0 = LayerModelEval(layers0,zd);
plot(t0,-zd,'-r'), hold on, plot(layers0(N+1:2*N),-layers0(1:N),'or'), plot(mask*conv(t0,kern),-zd,'--r');
plot(obs(:,2),-zd,'--b'), plot(trueProfile,-zd,'-b'), hold off;
legend('Inverse Model','','Blur of Model','Observed','True');
title(N+" Layers, Initial Model");
subplot(1,2,2);
t = LayerModelEval(layers,zd);
plot(t,-zd,'-r'), hold on, plot(layers(Nf+1:2*Nf),-layers(1:Nf),'or'), plot(mask*conv(t,kern),-zd,'--r');
plot(obs(:,2),-zd,'--b'), plot(trueProfile,-zd,'-b'), hold off;
legend('Inverse Model','','Blur of Model','Observed','True');
title(Nf+" Layers, Final Model");
