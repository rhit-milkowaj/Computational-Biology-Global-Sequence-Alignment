function F = alignmenting(F, gap_penalty, S, len1, len2)
%

for i=1:len1 % i is row index
    for j=1:len2 % j is column index
        down = F(i,j+1)+gap_penalty;
        across = F(i+1,j)+gap_penalty;
        diag = F(i,j) + S(i,j);
        
        F(i+1,j+1)=min([down,across,diag]);
        
    end
end



end