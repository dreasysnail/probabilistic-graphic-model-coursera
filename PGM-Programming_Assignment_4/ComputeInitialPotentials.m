%COMPUTEINITIALPOTENTIALS Sets up the cliques in the clique tree that is
%passed in as a parameter.
%
%   P = COMPUTEINITIALPOTENTIALS(C) Takes the clique tree skeleton C which is a
%   struct with three fields:
%   - nodes: cell array representing the cliques in the tree.
%   - edges: represents the adjacency matrix of the tree.
%   - factorList: represents the list of factors that were used to build
%   the tree. 
%   
%   It returns the standard form of a clique tree P that we will use through 
%   the rest of the assigment. P is struct with two fields:
%   - cliqueList: represents an array of cliques with appropriate factors 
%   from factorList assigned to each clique. Where the .val of each clique
%   is initialized to the initial potential of that clique.
%   - edges: represents the adjacency matrix of the tree. 
%
% Copyright (C) Daphne Koller, Stanford University, 2012


function P = ComputeInitialPotentials(C)

% number of cliques
N = length(C.nodes);

% initialize cluster potentials 
P.cliqueList = repmat(struct('var', [], 'card', [], 'val', []), N, 1);
P.edges = zeros(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%
% First, compute an assignment of factors from factorList to cliques. 
% Then use that assignment to initialize the cliques in cliqueList to 
% their initial potentials. 

% C.nodes is a list of cliques.
% So in your code, you should start with: P.cliqueList(i).var = C.nodes{i};
% Print out C to get a better understanding of its structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = length(C.factorList);
% set card 
allvars = [];
allcards = [];
for i = 1:F
    [allvars,~,ib] = union(allvars, C.factorList(i).var, 'stable');
    allcards = [allcards, C.factorList(i).card(ib)];
end  

for i = 1:N
    P.cliqueList(i).var = C.nodes{i};
    for v = P.cliqueList(i).var(:)'
        P.cliqueList(i).card = [P.cliqueList(i).card,allcards(allvars==v)];
    end
    P.cliqueList(i).val = ones(1,prod(P.cliqueList(i).card));
end
nextfactor = false;
for f = 1:F
    for i = 1:N
        if nextfactor
            nextfactor = false;
            if i~=1
                break
            end
        end
        if all(ismember(C.factorList(f).var, P.cliqueList(i).var))
            P.cliqueList(i) = FactorProduct(P.cliqueList(i),C.factorList(f));
            % assign to only one node
            nextfactor = true;
        end
    end
end

P.edges = C.edges;

end

