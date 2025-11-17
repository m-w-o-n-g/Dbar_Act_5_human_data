function [new_voltages]=Synthesize_Voltages(V,Jappl,T,K)
% This function synthesizes voltage patterns one would get from an o.n. set of CPs (such as trig)
% given measured voltages from a set of measured data arising from orthogonal patterns with the
% same number of linearly independent vectors (patterns) as that of the desired set

% Variables:
% L        number of electrodes
% K        number of CPs
% V        measured voltages (could be complex), row is electrode number, column is CP, size L x K
% Jappl       applied currents (could be complex), row is electrode number, column is CP, size L x K
% T        desired orthogonal current patterns (real)

% new_voltages        synthesized complex voltages from CPs T, size L x K


% row is electrode number
% column is CP
% Q(ell,k)=<T(:,k),Jappl(:,ell)>/<T(:,k),T(:,k)>, size = K x K


Q = zeros(K,K);
for ell=1:K
  for k=1:K
    Q(ell,k)=(T(:,k)'*Jappl(:,ell))/(T(:,k)'*T(:,k));
  end
end
S = inv(Q); 
new_voltages = V*S.';


