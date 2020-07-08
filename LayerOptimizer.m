function [layers,info] = LayerOptimizer(N,obs,kern,options)
% LayerOptimizer uses the Particle Swarm optimizer to fit a soil profile
% model to observed data via the pointspread function kern. Requires the
% image processing toolbox and the global optimization toolbox.
% 
% Jon Cooper
%
% Input:
%   N       - The number of initial layers in the model.
%   obs     - A Mx2 matrix with M soil resistance observations. The first
%               column is the (positive) depth vector, the second contains 
%               the observed resistances.
%   kern    - A vector representing a discrete pointspread function, used
%               in the context of 'blurry_data = conv(true_data,kern)'.
% Options : To be specified in name-value cell array,
%               e.g. options = {'AOITolType','or';'LOOpTol',5;'AOIpTol',5};
%   layers0 - A 2N column vector defining a soil layer model with N layers.
%               The first N elements are the layer positions, given in
%               absolute positive depth, the second N elements are the
%               respective layer resistances.
%   HybridFcn - Used in PSO options. ({'fmincon'}, [], PSO HybridFcn options)
%   LOOTolType : ('iterative','total','absolute','and',{'or'})
%           iterative - Exit LOO if misfit increases by more than LOOpTol
%                           percent compared to the previous iteration.
%           total - Exit LOO if misfit increases by more than LOOpTol
%                       percent compared to the initial misfit.
%           absolute - Exit LOO if misfit score is greater than LOOTol.
%           and - Exit LOO if the conditions for both 'total' and
%                   'absolute' are met.
%           or - Exit LOO if either of the conditions for 'total' and
%                   'absolute' are met.
%   LOOpTol - Percent allowed change in misfit for LOO, if used. {5}
%   LOOTol - Maximum absolute misfit score. {15}
%   AOITolType : ('iterative','absolute',{'and'},'or')
%           iterative - Only accept adding a new layer if it decreases the
%                           misfit by more than AOIpTol percent, compared
%                           to the previous iteration.
%           absolute - Only accept adding a new layer if the misfit is
%                           below AOITol.
%           and - Only accept a new layer if the conditions for both
%                   'iterative' and 'absolute' are met.
%           or - Only accept a new layer if either of the the conditions for
%                   'iterative' and 'absolute' are met.
%   AOIpTol - Percent allowed change in misfit for AOI, if used. Should be
%               strictly greater than LOOpTol, if both are used. {10}
%   AOITol - Minimum absolute misfit score. Should be strictly less than
%               LOOTolType, if both are used. {10}
%   MisfitNorm - The norm type to be used in calculating the misfit.
%                   ({2}, vector norm options)
%   ConvTol - The maximum allowed change in the layer model to stop the
%               main algorithm loop (in the inf-sense). {1e-3}
%   MIN - Smallest possible layer resistance. {0.1}
%   MAX - Largest possible layer resistance. {10}
%   Plt - Boolean for plotting misfits in AOI and LOO. {0}
% Output:
%   layers - a soil layer model (see layers0 description)

%% Define default options & read user-defined options
layers0 = [];
HybridFcn = 'fmincon';
LOOTolType = 'or';
LOOpTol = 5;
LOOTol = 15;
AOITolType = 'and';
AOIpTol = 10;
AOITol = 10;
MisfitNorm = 2;
ConvTol = 1e-3;
MIN = .1;
MAX = 10;
Plt = 0;

% rewrite default parameters if needed
if nargin == nargin(mfilename)
  for j = 1:size(options,1), eval([options{j,1},'= options{j,2};']); end
end

%% Unpack variables and initialize layer model

zd = obs(:,1); qcMeas = obs(:,2);
psf = kern(:);
LB = [zeros(N,1);zeros(N,1)+MIN];
UB = [max(zd)*ones(N,1);MAX*ones(N,1)];

% define the misfit function
mask = eye(size(obs,1));
mask = [zeros(size(obs,1),floor(length(kern)/2)),mask,zeros(size(obs,1),ceil(length(kern)/2)-1)];
Misfit=@(l) norm(mask*conv(LayerModelEval(l,zd),psf)-qcMeas,MisfitNorm);

% define the initial layer model, if not given, and set properties
nPop = min(100,10*2*N);
PSOoptions = optimoptions('particleswarm','SwarmSize',nPop,'HybridFcn',HybridFcn,'Display','off');
if isempty(layers0)
    layers0 = particleswarm(Misfit,2*N,LB,UB,PSOoptions)';
end
pos = layers0(1:N); res = layers0(N+1:end);
[pos,order] = sort(pos); res = res(order);
layers0 = [pos;res]; layers0(1) = 0;

%% Main Loop

layers = layers0;
layersPrev = layers + 1;
while length(layersPrev) ~= length(layers) || norm((layersPrev-layers)./UB,inf) > ConvTol
    layersPrev = layers;
    layers = LOO(layers,Misfit,LOOTol,LOOpTol,LOOTolType,Plt); Nc = length(layers)/2;
    
    PSOoptions.InitialSwarmMatrix = [layers';layers'];
    LB = [zeros(Nc,1);zeros(Nc,1)+MIN];
    UB = [max(zd)*ones(Nc,1);MAX*ones(Nc,1)]; % shrink if necessary
    
    layers = particleswarm(Misfit,2*Nc,LB,UB,PSOoptions)';
    pos = layers(1:Nc); res = layers(Nc+1:end);
    [pos,order] = sort(pos); res = res(order);
    layers = [pos;res]; layers(1) = 0;
    
    layers = AOI(layers,Misfit,AOITol,AOIpTol,AOITolType,UB,PSOoptions,Plt); Nc = length(layers)/2;
%     LB = [zeros(Nc,1);zeros(Nc,1)+MIN]; % not used
    UB = [max(zd)*ones(Nc,1);MAX*ones(Nc,1)]; % grow if necessary
end

%% Post processing

layers = LOO(layers,Misfit,LOOTol,LOOpTol,LOOTolType,Plt); Nc = length(layers)/2;
pos = layers(1:Nc); res = layers(Nc+1:end);
[pos,order] = sort(pos); res = res(order);
layers = [pos;res]; layers(1) = 0;

info = struct('Nfinal',Nc,'layers0',layers0,'mask',mask,'Misfit',Misfit,...
    'InitialMisfit',Misfit(layers0),'FinalMisfit',Misfit(layers));

end