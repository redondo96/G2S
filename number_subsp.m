function [nsub3, nsub2] = number_subsp(n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% We prioritize spaces of dimension 3. This could be changed

if mod(n,3)==0
    nsub3 = n/3;
    nsub2 = 0;
elseif mod(n,3)==1
    nsub3 = (n-4)/3;
    nsub2 = 2;
elseif mod(n,3)==2
    nsub3 = (n-2)/3;
    nsub2 = 1;
end

end

