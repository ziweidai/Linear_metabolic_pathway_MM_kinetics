function [J,e] = maxFluxEfficiency(Sin,Sout,kcat,Ks,Kp,Keq)
nrxn = length(kcat);
x0 = ones(nrxn,1)/nrxn;
A = ones(1,nrxn);
b = 1;
lb = zeros(nrxn,1);
ub = ones(nrxn,1);

%options = optimoptions(@fmincon,'MaxIterations',5000, ...
%    'OptimalityTolerance',1e-10);
e = fmincon(@(x)-SS_Linear(Sin,Sout,kcat.*x,Ks,Kp,Keq),x0,[],[],A,b,lb,ub);
J = SS_Linear(Sin,Sout,kcat.*e,Ks,Kp,Keq);
    
%{
options = optimoptions('fmincon',SpecifyConstraintGradient=true);
loge = fmincon(@(x)-SS_Linear(Sin,Sout,kcat.*(10.^x),Ks,Kp,Keq),log10(x0), ...
    [],[],[],[],-Inf*ones(nrxn,1),zeros(nrxn,1),@upper_bound_e_log);
J0 = SS_Linear(Sin,Sout,kcat.*e,Ks,Kp,Keq);
J1 = SS_Linear(Sin,Sout,kcat.*10.^loge,Ks,Kp,Keq);
[J0 J1]
if J0 > J1
    J = J0;
else
    J = J1;
    e = 10.^loge;
end
%}
end