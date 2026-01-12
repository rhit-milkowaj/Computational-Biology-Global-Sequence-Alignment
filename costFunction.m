function cost = costFunction(i, j, scoreMatrix)
% returns the cost for pairing the sequence elements i and j based on the
% scoring matrix

index_i = 0;
index_j = 0;

if i=='a'
    index_i = 1;
elseif i=='t'
    index_i = 2;
elseif i =='g'
    index_i = 3;
else
    index_i = 4;
end

if j=='a'
    index_j = 1;
elseif j=='t'
    index_j = 2;
elseif j =='g'
    index_j = 3;
else
    index_j = 4;
end

cost = scoreMatrix(index_i, index_j);
return




