function layers = LOO(layers,Misfit,TOL,pTOL,tolType,p)
% LOO (Leave-One-Out) repeatedly tries to remove the least-contributing
% layers until the defined tolerance for misfit increase is met. Does not
% re-optimize the layers.
%
% Jon Cooper
%
% Input:
%   layers  - layer model vector
%   Misfit  - misfit function
%   TOL     - either a percent or a fixed value, depends on tolType
%   pTOL    - percent allowed change in misfit, if used
%   tolType - 'iterative', 'total', 'absolute', 'and', or 'or'; see
%               LayerOptimizer options for more details
%   p       - plots misfits if greater than 2

base = Misfit(layers);
if strcmp(tolType,'iterative') || strcmp(tolType,'total')
    TOL = base*(pTOL/100+1);
elseif strcmp(tolType,'and')
    TOL = min([TOL,base*(pTOL/100+1)]);
elseif strcmp(tolType,'or')
    TOL = max([TOL,base*(pTOL/100+1)]);
end
N = length(layers)/2; Nc = N;% number of current layers
while 1 % loop condition not actually necessary; breaks inside
    if N == 1
        break
    end
    
    if strcmp(tolType,'iterative') % recompute the tolerance if desired
        TOL = Misfit(layers)*(pTOL/100+1);
    end
    
    misfits = zeros(Nc,1); % compute misfits for leave-1-out models
    for i = 1:Nc
        misfits(i) = Misfit(layers([1:i-1,i+1:Nc+i-1,Nc+i+1:end]));
    end
    if p > 2 && Nc == N % if first iteration, plot initial misfits
        figure, plot(sort(misfits)), hold on;
        yline(base,'--k'); yline(TOL,'--r');
        legend('Misfit','Original Misfit','Tolerance','location','NorthWest');
        title('Initial LOO Misfits');
    end
    
    if ~any(misfits <= TOL) % if jump in misfit is too great, exit
        if p > 2 % plot final misfits
            figure, plot(sort(misfits)), hold on;
            yline(base,'--k'); yline(TOL,'--r');
            legend('Misfit','Original Misfit','Tolerance','location','NorthWest');
            title('Final LOO Misfits');
        end
        break
    else
        [~,i] = min(misfits); % otherwise, delete the layer that has the least-worst impact on misfit when removed
        layers = layers([1:i-1,i+1:Nc+i-1,Nc+i+1:end]);
    end
    
    Nc = length(layers)/2; % number of layer left
    
    if Nc == 1
        break
    end
end
end