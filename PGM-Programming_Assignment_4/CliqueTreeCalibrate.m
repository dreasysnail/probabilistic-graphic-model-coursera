%CLIQUETREECALIBRATE Performs sum-product or max-product algorithm for 
%clique tree calibration.

%   P = CLIQUETREECALIBRATE(P, isMax) calibrates a given clique tree, P 
%   according to the value of isMax flag. If isMax is 1, it uses max-sum
%   message passing, otherwise uses sum-product. This function 
%   returns the clique tree where the .val for each clique in .cliqueList
%   is set to the final calibrated potentials.
%
% Copyright (C) Daphne Koller, Stanford University, 2012

function P = CliqueTreeCalibrate(P, isMax)


% Number of cliques in the tree.
N = length(P.cliqueList);

% Setting up the messages that will be passed.
% MESSAGES(i,j) represents the message going from clique i to clique j. 
MESSAGES = repmat(struct('var', [], 'card', [], 'val', []), N, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We have split the coding part for this function in two chunks with
% specific comments. This will make implementation much easier.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% YOUR CODE HERE
% While there are ready cliques to pass messages between, keep passing
% messages. Use GetNextCliques to find cliques to pass messages between.
% Once you have clique i that is ready to send message to clique
% j, compute the message and put it in MESSAGES(i,j).
% Remember that you only need an upward pass and a downward pass.
%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isMax
    for i=1:length(P.cliqueList)
        P.cliqueList(i).val = log(P.cliqueList(i).val);
    end
end

 while true
  [fromN,toN] = GetNextCliques(P,MESSAGES);  
   if fromN==0 || toN == 0
       break
   end
   idx = find( P.edges(:,fromN)==1 );
   idx = setdiff(idx,toN);
   % fromN -> toN = sum(phi_fromN * prod )
   message = P.cliqueList(fromN);
   % FACTORS.Joint = prod(factors)
   for i = 1:length(idx)
       if isMax
           message = FactorSum(message,MESSAGES(idx(i),fromN));
       else
           message = FactorProduct(message,MESSAGES(idx(i),fromN));
       end
   end
   if isMax
       message = FactorMaxMarginalization(message,setdiff(P.cliqueList(fromN).var,P.cliqueList(toN).var));
   else
       message = FactorMarginalization(message,setdiff(P.cliqueList(fromN).var,P.cliqueList(toN).var));
       message.val = message.val./sum(message.val);
   end
   MESSAGES(fromN,toN) = message;  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%
% Now the clique tree has been calibrated. 
% Compute the final potentials for the cliques and place them in P.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N
    if isMax
       P.cliqueList(i) = ComputeJointSumDistribution([P.cliqueList(i);MESSAGES(P.edges(:,i)==1,i)]);
    else
       P.cliqueList(i) = ComputeJointDistribution([P.cliqueList(i);MESSAGES(P.edges(:,i)==1,i)]);
    end    
end


return
