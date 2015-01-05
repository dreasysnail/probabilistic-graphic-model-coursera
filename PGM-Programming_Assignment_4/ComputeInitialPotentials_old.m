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
for i = 1:N
    P.cliqueList(i).var = C.nodes{i};
    % set card (factorlist is ordered)
    for v = P.cliqueList(i).var
        P.cliqueList(i).card = [P.cliqueList(i).card,C.factorList(v).card(1)];
    end
    P.cliqueList(i).val = ones(P.cliqueList(i).card);

%    P.cliqueList(i).val = P.cliqueList(i).val(:);
    
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
            %fill
            [difSet,idiff] = setdiff(P.cliqueList(i).var,C.factorList(f).var);
            C.factorList(f).var = [C.factorList(f).var,difSet];
            C.factorList(f).card = [C.factorList(f).card,P.cliqueList(i).card(idiff)];
            times = prod(P.cliqueList(i).card(idiff));
            C.factorList(f).val = repmat(C.factorList(f).val',times,1);
            %reorder
            [~,idx] = sort(C.factorList(f).var);
            tempMat = reshape(C.factorList(f).val,C.factorList(f).card);
            tempVal = permute(tempMat,idx);
            P.cliqueList(i).val = P.cliqueList(i).val .* tempVal;
            % assign to only one node
            nextfactor = true;
        end
    end
end

for i = 1:N
    P.cliqueList(i).val = P.cliqueList(i).val(:)';    
end
P.edges = C.edges;

end

