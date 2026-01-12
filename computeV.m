function [V] = computeV(s1_index,s2_index,string1,string2)
% s1_index = string 1 index
% s2_index = string 2 index

    if string1(s1_index) == string2(s2_index)
        V = -1;
    elseif (string1(s1_index) == 'a' && string2(s2_index) == 't') || ... 
           (string1(s1_index) == 't' && string2(s2_index) == 'a') || ...
           (string1(s1_index) == 'c' && string2(s2_index) == 'g') || ...
           (string1(s1_index) == 'g' && string2(s2_index) == 'c')
        V = 1;
    else
        V = 1;
    end

end

