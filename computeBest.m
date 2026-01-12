function [bestScore] = computeBest(string1,string2,gP)

sLength = length(string1);

F = zeros(sLength+1);

F(1,:) = 0:gP:(sLength*gP);
F(:,1) = 0:gP:(sLength*gP);


% Find matrix F
for i=2:sLength+1   % Loop through rows
    for j=2:sLength+1   % Loop through columns

        V = computeV(i-1,j-1,string1,string2);   % Compute V

        % Find F
        F1 = F(i-1,j)+gP;
        F2 = F(i,j-1)+gP;
        F3 = F(i-1,j-1)+V;

        F(i,j) = min([F1,F2,F3]);
    end

end

bestScore = F(sLength+1,sLength+1);   % Best F is the bottom right one
