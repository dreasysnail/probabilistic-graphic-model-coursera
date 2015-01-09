% function [nll, grad] = InstanceNegLogLikelihood(X, y, theta, modelParams)
% returns the negative log-likelihood and its gradient, given a CRF with parameters theta,
% on data (X, y). 
%
% Inputs:
% X            Data.                           (numCharacters x numImageFeatures matrix)
%              X(:,1) is all ones, i.e., it encodes the intercept/bias term.
% y            Data labels.                    (numCharacters x 1 vector)
% theta        CRF weights/parameters.         (numParams x 1 vector)
%              These are shared among the various singleton / pairwise features.
% modelParams  Struct with three fields:
%   .numHiddenStates     in our case, set to 26 (26 possible characters)
%   .numObservedStates   in our case, set to 2  (each pixel is either on or off)
%   .lambda              the regularization parameter lambda
%
% Outputs:
% nll          Negative log-likelihood of the data.    (scalar)
% grad         Gradient of nll with respect to theta   (numParams x 1 vector)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

function [nll, grad] = InstanceNegLogLikelihood(X, y, theta, modelParams)

    % featureSet is a struct with two fields:
    %    .numParams - the number of parameters in the CRF (this is not numImageFeatures
    %                 nor numFeatures, because of parameter sharing)
    %    .features  - an array comprising the features in the CRF.
    %
    % Each feature is a binary indicator variable, represented by a struct 
    % with three fields:
    %    .var          - a vector containing the variables in the scope of this feature
    %    .assignment   - the assignment that this indicator variable corresponds to
    %    .paramIdx     - the index in theta that this feature corresponds to
    %
    % For example, if we have:
    %   
    %   feature = struct('var', [2 3], 'assignment', [5 6], 'paramIdx', 8);
    %
    % then feature is an indicator function over X_2 and X_3, which takes on a value of 1
    % if X_2 = 5 and X_3 = 6 (which would be 'e' and 'f'), and 0 otherwise. 
    % Its contribution to the log-likelihood would be theta(8) if it's 1, and 0 otherwise.
    %
    % If you're interested in the implementation details of CRFs, 
    % feel free to read through GenerateAllFeatures.m and the functions it calls!
    % For the purposes of this assignment, though, you don't
    % have to understand how this code works. (It's complicated.)
    
    featureSet = GenerateAllFeatures(X, modelParams);

    % Use the featureSet to calculate nll and grad.
    % This is the main part of the assignment, and it is very tricky - be careful!
    % You might want to code up your own numerical gradient checker to make sure
    % your answers are correct.
    %
    % Hint: you can use CliqueTreeCalibrate to calculate logZ effectively. 
    %       We have halfway-modified CliqueTreeCalibrate; complete our implementation 
    %       if you want to use it to compute logZ.
    
    nll = 0;
    grad = zeros(size(theta));
    %%%
    % Your code here:
    % singleton: event(hiddenst, featureNum, obs)
    % factor:  when var == hiddenst, x, sum of feature
    [numVar, ~] = size(X);
    singletonFactors = repmat(struct('var', [],...
                              'card', modelParams.numHiddenStates,...
                              'val',  zeros(1,modelParams.numHiddenStates)), numVar, 1);
    pairFactors = repmat(struct('var', [],...
                              'card', modelParams.numHiddenStates*ones(1,2),...
                              'val',  zeros(1,modelParams.numHiddenStates^2)), numVar-1, 1);
    for v = 1:numVar
        singletonFactors(v).var = v;
    end
    c = containers.Map();
    makeKey = @(x) (num2str(x(1)*numVar+x(2)));
    for v = 1:numVar-1
        pairFactors(v).var = [v,v+1];
        c(makeKey(pairFactors(v).var)) = v;
    end
    
    for f = featureSet.features
        if length(f.var) == 1
            singletonFactors(f.var).val(f.assignment) = ...
                singletonFactors(f.var).val(f.assignment) + theta(f.paramIdx);   
        else 
            assert(length(f.var) == 2);
            idx = AssignmentToIndex(f.assignment, modelParams.numHiddenStates*ones(1,2));
            pairFactors(c(makeKey(f.var))).val(idx) = ...
                pairFactors(c(makeKey(f.var))).val(idx) + theta(f.paramIdx);
        end
        if isequal(y(f.var), f.assignment)
            nll = nll - theta(f.paramIdx);
        end
    end
    allFactors = [singletonFactors;pairFactors];
    for i = 1:length(allFactors)
        allFactors(i).val = exp(allFactors(i).val);
    end
    % 
    
    P = CreateCliqueTree(allFactors);
    [P, logZ] = CliqueTreeCalibrate(P, 0);
    nll = nll + logZ + modelParams.lambda / 2 * (theta * theta');
    
    grad = zeros(size(theta));
    %compute p(y'|x)  
    M = repmat(struct('var', [], 'card', [], 'val', []), numVar*2-1, 1);
    N = length(P.cliqueList);
    % singleton
    for i = 1:N
        for v = P.cliqueList(i).var(:)'
            if ~isempty(M(v).var)
                continue
            end
            M(v) = FactorMarginalization(P.cliqueList(i),setdiff(P.cliqueList(i).var,v));
            M(v).val = M(v).val/sum(M(v).val);
        end
    end
    % pairwise
    for i = 1:N
        M(i+numVar) = P.cliqueList(i);
        M(i+numVar).val = M(i+numVar).val/sum(M(i+numVar).val);
    end
    
    for f = featureSet.features
        %feature count
        grad(f.paramIdx) = grad(f.paramIdx) - isequal(f.assignment, y(f.var));
        %expect count
        if length(f.var) == 1
            expectCount = M(f.var).val(f.assignment);
        else
            idx = AssignmentToIndex(f.assignment, modelParams.numHiddenStates*ones(1,2));
            expectCount = M(f.var(1)+numVar).val(idx);
        end
        grad(f.paramIdx) = grad(f.paramIdx) + expectCount;
    end
    grad = grad + modelParams.lambda * theta;
    
    
    
end
