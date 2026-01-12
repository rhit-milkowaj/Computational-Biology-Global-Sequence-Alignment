
% dynamic programming scratch work


%% Try to replicate intro example


seq1 = 'attcga';
seq2 = 'aattga';

scoreMatrix = [0 1 4 4; % A T G C
               1 0 4 4;
               4 4 0 1;
               4 4 1 0];
gap_penalty = 2;

F = zeros([length(seq1)+1,length(seq2)+1]);
F(1,:) = 0:gap_penalty:length(seq1)*gap_penalty;
F(:,1) = 0:gap_penalty:length(seq2)*gap_penalty;
traceBack = [];

for i=2:length(seq1)+1 % i is row index
    for j=2:length(seq2)+1 % j is column index
        down = F(i-1,j)+gap_penalty;
        across = F(i,j-1)+gap_penalty;
        diagonal = F(i-1,j-1) + costFunction(seq1(j-1),seq2(i-1),scoreMatrix);

        [F(i,j),index]=min([down,across,diagonal]);
        if index == 1
            traceBack(i,j,:)=[i-1,j, index];
        elseif index == 2
            traceBack(i,j,:)=[i,j-1, index];
        else 
            traceBack(i,j,:)=[i-1,j-1, index];
        end
    end
end

traceBack(1,:,:) = [0,ones([1,length(seq1)]);0:length(seq1);2*ones([1,length(seq1)+1])]';
traceBack(:,1,:) = [0:length(seq1);0,ones([1,length(seq1)]);ones([1,length(seq1)+1])]';

revSeq1 = '';
revSeq2 = '';
while i>1 && j>1
    [i,j,revSeq1,revSeq2] = tracer(i,j,traceBack,seq1,seq2,revSeq1,revSeq2);
    
end

newSeq1 = reverse(revSeq1)
newSeq2 = reverse(revSeq2)
F(7,7)







%% Protein Alignment
close all;

gap_penalty = 0.3;

load proteins.mat

x_length = 44;
y_length = 44;
z_data = [];

for k=1:x_length
    for i=1:y_length
        protien_k = proteins(k);
        protien_k = protien_k{1};
        protien_i = proteins(i);
        protien_i = protien_i{1};

        S = get_Simularity_values(protien_k,protien_i);

        F = zeros([length(protien_k)+1,length(protien_i)+1]);
        F(:,1) = 0:gap_penalty:length(protien_k)*gap_penalty;
        F(1,:) = 0:gap_penalty:length(protien_i)*gap_penalty;
        

        F = alignmenting(F, gap_penalty, S, length(protien_k), length(protien_i));
        [y,x]=size(F);
        z_data(k,i) = F(y,x);

    end
end
figure(1)
pcolor(z_data)

%% finding best k for clustering
for ii=1:40

n = 44;
k = ii;
D = z_data;

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
x_square = reshape(gmSoln.x(1:n*n),n,n);

I = find(diag(x_square));
new_order = [];
Groups = {};

for i=1:k
    new_order = cat(2,new_order,find(x_square(I(i),:)));
    Groups(i) = {find(x_square(I(i),:))};
end

new_z_data = zeros(n);

for i=1:n
    new_z_data(i,:) = z_data(new_order(i),new_order);
end

figure(ii)
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


%% getting figure for specific k values

n = 44;
k = 10;
D = z_data;

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
x_square = reshape(gmSoln.x(1:n*n),n,n);

I = find(diag(x_square));
new_order = [];
Groups = {};

for i=1:k
    new_order = cat(2,new_order,find(x_square(I(i),:)));
    Groups(i) = {find(x_square(I(i),:))};
end

new_z_data = zeros(n);

for i=1:n
    new_z_data(i,:) = z_data(new_order(i),new_order);
end

figure(2)
pcolor(new_z_data)
Groups;

clusterScore = 0;
for i=1:k
    group_i = Groups(i);
    group_i = group_i{1};
    group_not_i = 1:44;
    group_not_i = group_not_i(~ismember(group_not_i,group_i));

    center = I(i);
    gi_scores = z_data(i,group_i);
    gni_scores = z_data(i,group_not_i);

    clusterScore = clusterScore + mean(gi_scores)/mean(gni_scores)/sqrt(k);


end
k
clusterScore



