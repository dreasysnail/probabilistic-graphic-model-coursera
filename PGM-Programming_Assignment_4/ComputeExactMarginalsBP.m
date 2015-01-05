%COMPUTEEXACTMARGINALSBP Runs exact inference and returns the marginals
%over all the variables (if isMax == 0) or the max-marginals (if isMax == 1). 
%
%   M = COMPUTEEXACTMARGINALSBP(F, E, isMax) takes a list of factors F,
%   evidence E, and a flag isMax, runs exact inference and returns the
%   final marginals for the variables in the network. If isMax is 1, then
%   it runs exact MAP inference, otherwise exact inference (sum-prod).
%   It returns an array of size equal to the number of variables in the 
%   network where M(i) represents the ith variable and M(i).val represents 
%   the marginals of the ith variable. 
%
% Copyright (C) Daphne Koller, Stanford University, 2012


function M = ComputeExactMarginalsBP(F, E, isMax)

% initialization
% you should set it to the correct value in your code
% M = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%
% Implement Exact and MAP Inference.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allvars = [];
for i = 1:length(F)
    allvars = union(allvars, F(i).var);
end   
M = repmat(struct('var', [], 'card', [], 'val', []), length(allvars), 1);
P = CreateCliqueTree(F, E);
P = CliqueTreeCalibrate(P,isMax);


N = length(P.cliqueList);
% not so efficient: should pick the optimal node; should build a map
for i=1:N
    for v = P.cliqueList(i).var(:)'
        if ~isempty(M(v).var)
            continue
        end
        if isMax
            M(v) = FactorMaxMarginalization(P.cliqueList(i),setdiff(P.cliqueList(i).var,v));
        else
            M(v) = FactorMarginalization(P.cliqueList(i),setdiff(P.cliqueList(i).var,v));
            M(v).val = M(v).val/sum(M(v).val);
        end
    end
end



end
