function out=simp1D(func,a,b,NX)

% This function will integrate the function 'func' from a to b with NX
% number of intervals

% Make sure NX is even
NX=2*ceil(NX/2);

h=(b-a)/NX;
xg=a:h:b;
U=feval(func,xg);
ixo = 2:2:NX;
ixe = 3:2:NX-1;

out=(h/3)*(U(1)+2*sum(U(ixe))+4*sum(U(ixo))+U(NX+1));