function layers = LOOlog(layers,Misfit,lTOL,rTOL)
% LOOLog (Leave-One-Out) is quite different from the original LOO
% algorithm. Now, we look at all pairs of adjacent layers, check if they
% are too close in either location or resistance, and remove the one that
% contributes least to the misfit.
%
% Jon Cooper
%
% Input:
%   layers  - layer model vector
%   Misfit  - misfit function
%   lTOL    - minimum layer thickness
%   rTOL    - minimum difference in layer resistance

N = length(layers)/2; Nc = N;% number of current layers
if N == 1
    return
end

j = 0;
for i = 1:Nc-1
    i = i-j; % intentional. j keeps track of layers removed. OK since i resets
    if layers(i+1)-layers(i) < lTOL || abs(layers(Nc+i+1)-layers(Nc+i)) < rTOL % check if violated
        misfitA = Misfit(layers([1:i-1,i+1:Nc+i-1,Nc+i+1:end]));
        misfitB = Misfit(layers([1:i,i+2:Nc+i,Nc+i+2:end]));
        
        % remove least-contributing layer of the two
        if misfitA < misfitB
            layers = layers([1:i-1,i+1:Nc+i-1,Nc+i+1:end]);
        else
            layers = layers([1:i,i+2:Nc+i,Nc+i+2:end]);
        end
        Nc = Nc-1; j = j+1;
    end
end
end