%% BMTH413 Project 2 - Dynamic Programming

% Written by Taylor Donen & Angela Milkowski


%% Gurobi - kMeans calculation
addpath('C:\gurobi1300\win64\matlab\');

dnaData = importdata('dnadata_project2.xlsx');

gP = 1;   % gP = Gap penalty = can adjust 
n = length(dnaData);   % Number of sequences
scoreMatrix = zeros(n);

for i=1:n   % Outer loop: rows 
    for j=1:n  % Inner loop: columns
        
        string1 = char(dnaData(i));
        string2 = char(dnaData(j));

        bestScore = computeBest(string1,string2,gP);
        scoreMatrix(i,j) = bestScore;

    end
end

% k-means using Gurobi
k = 15;
D = scoreMatrix;

A = zeros(n+1,n*n);
B = repmat(-eye(n),1,n);

for i=1:n
   A(i,n*(i-1)+1:n*i) = ones(1,n);
   A(n+1,n*(i-1)+i) = 1;
   B(i,n*(i-1)+i) = n-k+1;
end

% Gurobi requires standard form :-(
A = [A zeros(n+1,n); B -eye(n)];

% gm is a Gurobi Model, define Gurobi model
gm.modelsense = 'min';
gm.obj = [D(:); zeros(n,1)];
gm.A = sparse(A);
gm.rhs = [ones(n,1);k;zeros(n,1)];
gm.vtype = [repmat('B',n*n,1);repmat('C',n,1)];
gm.sense = '=';

% Solve
gmSoln = gurobi(gm);
x = reshape(gmSoln.x(1:n*n),n,n);

%% Graphing Heatmap

nonzeroInd = find(diag(x)==1);
numGroups = length(nonzeroInd);

groups = cell(numGroups,1);

for i=1:numGroups
    groups{i} = find(x(nonzeroInd(i), :) == 1);
end

plotDP = [];

for i=1:numGroups
    row = [];
    for j=1:numGroups
        row = [row, D(groups{i}, groups{j})];
    end

    plotDP = [plotDP; row];
end
          
% figure(1);
% imagesc(plotDP)      % Scale data and display as image
% colorbar                  % Show scale of scores
% 
% xlabel('Sequence index j')
% ylabel('Sequence index i')
% title('Map of Best Alignment Scores (lower = better)')

%% Find best k

for ii=1:40

n = length(dnaData);   % Number of sequences
k = ii;
D = scoreMatrix;
z_data = D;

A = zeros(n+1,n*n);
B = repmat(-eye(n),1,n);
for i=1:n
   A(i,n*(i-1)+1:n*i) = ones(1,n);
   A(n+1,n*(i-1)+i) = 1;
   B(i,n*(i-1)+i) = n-k+1;
end
% Gurobi requires standard form :-(
A = [A zeros(n+1,n); B -eye(n)];

% gm is a Gurobi Model
gm.modelsense = 'min';
gm.obj = [D(:); zeros(n,1)];
gm.A = sparse(A);
gm.rhs = [ones(n,1);k;zeros(n,1)];
gm.vtype = [repmat('B',n*n,1);repmat('C',n,1)];
gm.sense = '=';

gmSoln = gurobi(gm);
x = reshape(gmSoln.x(1:n*n),n,n);

I = find(diag(x));
new_order = [];
Groups = {};

for i=1:k
    new_order = cat(2,new_order,find(x(I(i),:)));
    Groups(i) = {find(x(I(i),:))};
end

new_z_data = zeros(n);

for i=1:n
    new_z_data(i,:) = z_data(new_order(i),new_order);
end

figure(2)
pcolor(new_z_data)
Groups;
lengths = [];

clusterScore = 0;
for i=1:k
    group_i = Groups(i);
    group_i = group_i{1};
    group_not_i = 1:n;
    group_not_i = group_not_i(~ismember(group_not_i,group_i));

    center = I(i);
    gi_scores = z_data(i,group_i);
    gni_scores = z_data(i,group_not_i);

    clusterScore = clusterScore + (length(group_i)-n/k)^3 + ...
                                (median(gi_scores)/median(gni_scores))^4 * ...
                                k;

end
k;
clusterScore;

ClusterScores(k) = clusterScore



end

[min_score,best_k] = min(ClusterScores)