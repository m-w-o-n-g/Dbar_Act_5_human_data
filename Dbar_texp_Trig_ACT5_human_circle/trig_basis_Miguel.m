function [E,e,ne] = trig_basis_Miguel(L)
% Computes the trig basis
% [E,e]=trig_basis(# of electrodes)
% E is a matrix of sines and cosines
% e is E normalized w.r.t. ||.||2

dtheta=2*pi/L;
% dx=r*dtheta;
theta_l=[dtheta: dtheta: dtheta*L].';

for n=1:L/2
	E(:,n)=cos(n*theta_l);
end

for n=1:L/2-1
	E(:,n+L/2)=sin(n*theta_l);
end
% ne=sqrt(diag(dx*E.'*E));
 ne=sqrt(diag(E.'*E));
for n=1:L-1
    e(:,n)=E(:,n)./ne(n);
end
