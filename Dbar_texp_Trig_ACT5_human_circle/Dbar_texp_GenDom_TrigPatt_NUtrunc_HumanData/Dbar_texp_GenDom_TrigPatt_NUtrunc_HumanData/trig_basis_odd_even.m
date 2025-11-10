function [E,e,ne] = trig_basis_odd_even(L,Ldiv2)
% Computes the trig basis
% [E,e]=trig_basis(# of electrodes)
% E is a matrix of sines and cosines
% e is E normalized w.r.t. ||.||2

% This works if the number L of electrodes is odd or even
% Ldiv2 is either floor(L/2) or ceil(L/2)

dtheta=2*pi/L;

theta_l=[dtheta: dtheta: dtheta*L].';

if (Ldiv2==ceil(L/2))  % in this case we will have one more cosine pattern than sines
E = zeros(L,L-1);

for n=1:Ldiv2
	E(:,n)=cos(n*theta_l);
end

for n=1:Ldiv2-1
	E(:,n+Ldiv2)=sin(n*theta_l);
end
               
else % In the case of floor, there will be one more sine than cosines
% Note that if L is even it doesn't matter, there will the same number of each
for n=1:Ldiv2
   E(:,n)=cos(n*theta_l);
end
               
for n=1:Ldiv2+1
   E(:,n+Ldiv2)=sin(n*theta_l);
end

end % end of if statement

 ne=sqrt(diag(E.'*E));
              
for n=1:L-1
    e(:,n)=E(:,n)./ne(n);
end
