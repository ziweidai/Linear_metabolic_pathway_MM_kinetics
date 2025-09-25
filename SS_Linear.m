function [J,met_conc] = SS_Linear(Sin,Sout,kcat,Ks,Kp,Keq)

%Calculate steady state metabolite concentrations based on input/output
%substrate concentrations and enzyme parameters
%Calculate flux control coefficients (FCCs), flux (J) and reaction free energies
%based on the steady state metabolite concentrations
%-----------------------------------------------------------------------
%Calculate steady state concentrations
n=length(kcat)-1;
%{
y0=((Sin*Sout)^0.5)*ones(n,1);
options = optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-10);
%y=fsolve(@(x) dxdt_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq),y0,options);
[y,resnorm] = lsqnonlin(@(x) dxdt_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq),y0,zeros(n,1),[],options);
%}
y0=log10(((Sin*Sout)^0.5)*ones(n,1));
options = optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-10);
%y=fsolve(@(x) dxdt_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq),y0,options);
[y,resnorm] = lsqnonlin(@(x) dxdt_Linear(10.^x,Sin,Sout,kcat,Ks,Kp,Keq), ...
    y0,[],[],options);
y = 10.^y;
met_conc = y;
if resnorm>1e-6  
    %resnorm
    J=-1000;
else    
    %Calculate flux
    JVec=zeros(n+1,1);
    S_ext=[Sin;y(:);Sout];
    for i=1:n+1
        JVec(i)=MM(S_ext(i),S_ext(i+1),kcat(i),Ks(i),Kp(i),Keq(i));
    end
    J=mean(JVec);
    if (max(JVec)-min(JVec))/J>0.01
        %JVec
        %max(JVec)-min(JVec)
        J=-1000;
    end
end
