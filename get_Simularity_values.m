function S = get_Simularity_values(D1,D2)
%imports the diagonal matricies containing the eigenvalues and returns the
%simularity scores matrix

len1 = length(D1);
len2 = length(D2);
S=[];

for i=1:len1
    for j=1:len2
        S(i,j) = abs(D1(i)-D2(j))/abs(D1(i)+D2(j));
    end
end


end