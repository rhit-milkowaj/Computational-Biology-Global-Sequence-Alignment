function [i,j,rev1,rev2] = tracer(i,j,traceBack, seq1,seq2,rev1,rev2)
%
index = traceBack(i,j,3);
if index == 2
    rev2=insertAfter(rev2,length(rev2),'_');
    rev1(end+1)=seq1(j-1); % this part is right I think
elseif index == 1
    rev2(end+1)=seq2(i-1); % this part is right I think
    rev1=insertAfter(rev1,length(rev1),'_');
else
    rev1(end+1)=seq1(j-1);
    rev2(end+1)=seq2(i-1);
end
temp_i = i;
i=traceBack(i,j,1);
j=traceBack(temp_i,j,2);

end