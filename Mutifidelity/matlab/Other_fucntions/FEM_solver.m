function u = FEM_solver(grid_level,xi,eta)
%function u = FEM_solver(n,xi,eta)
% grid_level: integer, with # of grids: N = 2*10*2^level+1;
% xi: randno, scalar
% eta: randno, scalar
%

n = 10*2^grid_level;
N = 2*n+1;
x = linspace(0,2,N);
h = 2/(N-1);
a1 = 1*(1+xi);
a2 = 10*(1+eta);
f1 = @(x) -2-3*x^2;%zeros(N,1);
f2 = @(x) 1-6*x;%zeros(N,1);
% f(N)=1;
A = sparse(N,N);
b = sparse(N,1);
freenodes = setdiff(1:N,[1,N]);
for i = 1:n
    elem = [i,i+1];
    A(elem,elem) = A(elem,elem) + a1/h*[1,-1;-1,1];
    b(elem) = b(elem) + f1(sum(x(elem))/2)*h/2;
end

for i = n+1:N-1
    elem = [i,i+1];
    A(elem,elem) = A(elem,elem) + a2/h*[1,-1;-1,1];
    b(elem) = b(elem) + f2(sum(x(elem))/2)*h/2;
end

u = sparse(N,1);
u(1) = 0;
u(end) = 1;
b = b - A * u;
u(freenodes) = A(freenodes,freenodes)\b(freenodes);
end