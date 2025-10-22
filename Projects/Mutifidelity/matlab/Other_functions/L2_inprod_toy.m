function result = L2_inprod_toy(u1,u2,h)
%function result = L2_inprod_toy(u1,u2,h)--inner product of u1 & u2
%INPUT: 
% u1, u2: vectors of function value
% h: mesh size
%
%OUTPUT: 
% result:scalar. L2 inner product of u1 & u2.
%
% Last modified: Feb-14, 2025

result = sum(h*u1.*u2);
end