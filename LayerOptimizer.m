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
%   LOOTol - Maximum absolute misfit score. {1}
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
%               LOOTol, if both are used. {1}
%   MisfitNorm - The norm type to be used in calculating the misfit.
%                   ({2}, vector norm options)
%   ConvTol - The maximum allowed change in the layer model to stop the
%               main algorithm loop (in the inf-sense). {1e-3}
%   MIN - Smallest possible layer resistance. {0.01}
%   MAX - Largest possible layer resistance. {10}
%   Plt - Integer for controlling what is plotted. 0 for none, 1 for initial,
%           coordinate descent, and final, 2 for LOO, AOI, and particleswarm
%           outputs within each iteration, 3 for all (including LOO & AOI misfits). {0}
%   MaxIter - Maximum number of iterations for the main loop. {10}
%   CD - Boolean for using coordinate descent on the initial layer model. {1}
%   CDMaxIter - Maximum number iterations for coordinate descent. {20}
% Output:
%   layers - a soil layer model (see layers0 description)

%% Define default options & read user-defined options
layers0 = [];
HybridFcn = 'fmincon';
LOOTolType = 'iterative';
LOOpTol = 5;
LOOTol = 1;
AOITolType = 'and';
AOIpTol = 10;
AOITol = 1;
MisfitNorm = 2;
ConvTol = 1e-3;
MIN = .01;
MAX = 10;
Plt = 0;
MaxIter = 10;
CD = 1;
CDMaxIter = 20;

% rewrite default parameters if needed
if nargin == nargin(mfilename)
  for j = 1:size(options,1), eval([options{j,1},'= options{j,2};']); end
end

if ~isempty(layers0), N = length(layers0)/2; end

%% Unpack variables and initialize layer model

M = size(obs,1);
zd = obs(:,1); qcMeas = obs(:,2);
mzd = max(zd); mqc = MAX;
zd = zd/mzd; qcMeas = qcMeas/mqc;
psf = kern(:);
LB = [zeros(N,1)+min(zd);zeros(N,1)+MIN];
UB = ones(2*N,1);

% define the misfit function
mask = eye(M); n = 0; z = 0;
mask = [zeros(M,floor(length(kern)/2)),mask,zeros(M,ceil(length(kern)/2)-1)];
mask = diag([zeros(1,z),linspace(0,1,n),ones(1,M-n-z)].*[ones(1,M-n-z),linspace(1,0,n),zeros(1,z)])*mask;
Misfit=@(l) norm(mask*conv(LayerModelEval(l,zd),psf)-qcMeas,MisfitNorm);

% define the initial layer model, if not given, and set properties
nPop = min(100,10*2*N);
PSOoptions = optimoptions('particleswarm','SwarmSize',nPop,'HybridFcn',HybridFcn,'Display','off');
if isempty(layers0)
    layers0 = particleswarm(Misfit,2*N,LB,UB,PSOoptions)';
    pos = layers0(1:N); res = layers0(N+1:end);
    [pos,order] = sort(pos); res = res(order);
    layers0 = [pos;res]; layers0(1) = 0;
else
    pos = layers0(1:N); res = layers0(N+1:end);
    [pos,order] = sort(pos); res = res(order);
    layers0 = [min(zd)+pos/mzd;res/mqc]; layers0(1) = 0;
end
if Plt > 0, layerPlot(zd,psf,mask,layers0,qcMeas,"Initial"); end

% coordinate descent for initial model
if CD
    CDoptions = optimoptions('fminunc','Display','off','Algorithm','quasi-newton',...
        'HessUpdate','steepdesc','MaxFunctionEvaluations',100); % backdoor gradient descent
    for j = 1:CDMaxIter
        layers = layers0;
        for i = length(layers0):-1:2
            layers0(i) = fminunc(@(l) Misfit([layers0(1:(i-1));l;layers0((i+1):(2*N))]),layers0(i),CDoptions);
        end
        if norm(layers-layers0) < 1e-2, break, end
    end
    
    % error checking / correcting
    pos = layers0(1:N); res = layers0(N+1:end);
    pos(pos > max(zd)) = NaN;
    res = res(~isnan(pos)); pos = pos(~isnan(pos));
    res(res < MIN) = MIN; res(res > MAX) = MAX;
    layers0 = [pos;res];

    if Plt > 0, layerPlot(zd,psf,mask,layers0,qcMeas,"Coordinate Descent"); end
end


%% Main Loop

iter = 0;
layers = layers0;
layersPrev = layers + 1;
while iter < MaxIter && (length(layersPrev) ~= length(layers) || norm(layersPrev-layers,inf) > ConvTol)
    layersPrev = layers; Nc = length(layers)/2;
    
    PSOoptions.InitialSwarmMatrix = [layers';layers'];
    layers = particleswarm(Misfit,2*Nc,LB,UB,PSOoptions)';
    if Plt > 1, layerPlot(zd,psf,mask,layers,qcMeas,"PSO A"); end
    
    layers = LOO(layers,Misfit,LOOTol,LOOpTol,LOOTolType,Plt); Nc = length(layers)/2;
    if Plt > 1, layerPlot(zd,psf,mask,layers,qcMeas,"LOO"); end
    
    PSOoptions.InitialSwarmMatrix = [layers';layers'];
    LB = [zeros(Nc,1)+min(zd);zeros(Nc,1)+MIN];
    UB = ones(2*Nc,1); % shrink if necessary
    
    layers = particleswarm(Misfit,2*Nc,LB,UB,PSOoptions)';
    pos = layers(1:Nc); res = layers(Nc+1:end);
    [pos,order] = sort(pos); res = res(order);
    layers = [pos;res]; layers(1) = 0;
    if Plt > 1, layerPlot(zd,psf,mask,layers,qcMeas,"PSO B"); end
    
    layers = AOI(layers,Misfit,AOITol,AOIpTol,AOITolType,UB,PSOoptions,Plt);
    if Plt > 1, layerPlot(zd,psf,mask,layers,qcMeas,"AOI"); end
    
    iter = iter+1;
end

%% Warnings
errorfunc = mask*conv(LayerModelEval(layers,zd),psf)-qcMeas;
error = norm(errorfunc,MisfitNorm);
if error > .2
    warning("High misfit: layer(s) likely averaged");
    
    % moving average, local L2 approximations
    zs = 1:length(zd);
    pos = (errorfunc(2:end) < 0) & (errorfunc(1:end-1) >= 0);
    neg = (errorfunc(2:end) > 0) & (errorfunc(1:end-1) <= 0);
    zs = zs(pos | neg); zs = [1;zs';length(zd)];
    
    for i = 1:length(zs)-2
        start = zs(i); fin = zs(i+2); width = zd(fin)-zd(start);
        eval1 = abs(sum(errorfunc(start:fin)))/width;
        eval2 = sum(abs(errorfunc(start:fin)))/width;

        if eval2/eval1 > 2 && eval2 > 10
            warning("Missing/incorrect layers likely between zd="+...
                string(zd(start)*mzd)+" and zd="+string(zd(fin)*mzd));
        end
    end
end


%% Post processing

layers = LOO(layers,Misfit,LOOTol,LOOpTol,LOOTolType,Plt); Nc = length(layers)/2;
pos = layers(1:Nc); res = layers(Nc+1:end);
[pos,order] = sort(pos); res = res(order);
layers = [pos;res]; layers(1) = 0;
fmisfit = Misfit(layers); imisfit = Misfit(layers0);
if Plt > 0, layerPlot(zd,psf,mask,layers,qcMeas,"Final"); end
layers = [(pos-min(zd))*mzd;res*mqc]; layers(1) = 0;
pos = layers0(1:N); res = layers0(N+1:end);
layers0 = [(pos-min(zd))*mzd;res*mqc]; layers0(1) = 0;

info = struct('Nfinal',Nc,'layers0',layers0,'mask',mask,'Misfit',Misfit,...
    'InitialMisfit',imisfit,'FinalMisfit',fmisfit);



end

function layerPlot(zd,psf,mask,model,meas,str)

figure
N = length(model)/2;
t = LayerModelEval(model,zd);
plot(t,-zd,'-r'), hold on, plot(model(N+1:2*N),-model(1:N),'or')
plot(mask*conv(t,psf),-zd,'--r'), plot(meas,-zd,'--b');
legend('Inverse Model','','Blur of Model','Observed');
title(str+", "+N+" Layers");

end