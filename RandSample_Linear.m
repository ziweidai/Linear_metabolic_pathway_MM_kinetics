nrxn=10;
nmet=nrxn-1;
nsamples=10000;
Keq=10.^normrnd(1,1,nrxn,nsamples);
kcat=10.^normrnd(1,1,nrxn,nsamples);
Ks=10.^normrnd(1,1,nrxn,nsamples);
Kp=10.^normrnd(1,1,nrxn,nsamples);
Sin=10.^normrnd(0,1,1,nsamples);
Sout=10.^normrnd(0,1,1,nsamples);

%count=0;
c = zeros(nsamples,1);
cadj = c;
J = zeros(nsamples,1);
eMat = zeros(nsamples,nrxn);
phiMat = zeros(nsamples,nrxn-1);
ae2rMat = phiMat;
saturationMat = eMat;
for i=1:nsamples
    i
    dG_overall=log(Sout(i)/Sin(i)/prod(Keq(:,i)));
    if dG_overall>0
        b=Sin(i);
        Sin(i)=Sout(i);
        Sout(i)=b;
        kcat(:,i)=kcat(nrxn:-1:1,i);
        Ks(:,i)=Ks(nrxn:-1:1,i);
        Kp(:,i)=Kp(nrxn:-1:1,i);
        Keq(:,i)=Keq(nrxn:-1:1,i);
    end
    [J(i),e] = maxFluxEfficiency(Sin(i),Sout(i),kcat(:,i),Ks(:,i),Kp(:,i),Keq(:,i));    
    [c(i),cadj(i),sat,phi,ae2r] = corr_Keq_aEratio(Keq(:,i),kcat(:,i),Ks(:,i),Kp(:,i),Sin(i),Sout(i),e);
    eMat(i,:) = reshape(e,1,nrxn);
    phiMat(i,:) = reshape(phi,1,nrxn-1);
    ae2rMat(i,:) = reshape(ae2r,1,nrxn-1);
    saturationMat(i,:) = reshape(sat,1,nrxn);
end 
figure;
hist(c(J>0),20);
xlabel("Spearman's \rho");
ylabel("Number of models");
title("Spearman correlation between K_i and a_i[E_i]^2/a_{i+1}[E_{i+1}]^2");

%% Analyze the relationship between the K-aE2 correlation and saturation terms
figure;
dscatter(mean(saturationMat(J>0,:),2),c(J>0));
xlabel("Mean saturation term, Michaelis-Menten kinetics");
ylabel("Spearman's \rho between K_i and a_i[E_i]^2/a_{i+1}[E_{i+1}]^2");
box on;
title("Relationship between saturation term and the K-aE^2 coupling");
hold on;
plot([0 1],[0 0],":");

%% Compare the spearman rho, unadjusted vs adjusted
figure;
subplot(1,2,1);
dscatter(c(J>0),cadj(J>0));hold on;plot([-1 1],[-1 1]);
xlabel("K_i vs a_i[E_i]^2/a_{i+1}[E_{i+1}]^2");
ylabel("K_i\Phi_i vs a_i[E_i]^2/a_{i+1}[E_{i+1}]^2");
box on;
title("Spearman's \rho between thermodynamic and kinetic terms")

subplot(1,2,2);
violinplot([c(J>0) cadj(J>0)],["K_i","K_i\Phi_i"], ...
    'ViolinColor',[87 127 219]/255);
box on;
title("Spearman's \rho with a_i[E_i]^2/a_{i+1}[E_{i+1}]^2");

%% Random simulation of relationship between eta=J/e and metabolite load
Sin = 100;
Sout = 1;
nrxn = 10;
nsamples = 10000;
Keq=10.^normrnd(1,1,nrxn,1);
Ks=10.^normrnd(1,1,nrxn,1);
Kp=10.^normrnd(1,1,nrxn,1);
kcat = ones(nrxn,1);
e = lhsdesign(nsamples,nrxn);
e = e./sum(e,2);
eta = zeros(nsamples,1);
sum_met_conc = eta;
for i = 1:nsamples
    [J(i),met_conc] = SS_Linear(Sin,Sout,kcat.*e(i,:)',Ks,Kp,Keq);
    sum_met_conc(i) = sum(met_conc);
end
[J_opt,e_opt] = maxFluxEfficiency(Sin,Sout,kcat,Ks,Kp,Keq);
[~,met_conc_opt] = SS_Linear(Sin,Sout,kcat.*e_opt,Ks,Kp,Keq);

% Plot relationship between eta and metabolite concentration
figure;
%subplot(1,2,1);
scatter(log(sum_met_conc(J>0)),log(J(J>0)),20,'Marker','+');
hold on;
scatter(log(sum(met_conc_opt)),log(J_opt),100,'filled');
box on;
fontsize(40,"points");
xlabel("log(total metabolite concentration)");
ylabel("log(flux efficiency)");
%title("No C-limitation");

Sin = 10;
for i = 1:nsamples
    [J(i),met_conc] = SS_Linear(Sin,Sout,kcat.*e(i,:)',Ks,Kp,Keq);
    sum_met_conc(i) = sum(met_conc);
end
[J_opt,e_opt] = maxFluxEfficiency(Sin,Sout,kcat,Ks,Kp,Keq);
[~,met_conc_opt] = SS_Linear(Sin,Sout,kcat.*e_opt,Ks,Kp,Keq);
%subplot(1,2,2);
scatter(log(sum_met_conc(J>0)),log(J(J>0)),20,'Marker','+');
hold on;
scatter(log(sum(met_conc_opt)),log(J_opt),100,'filled');
box on;
legend(["No C-lim random","No C-lim optimal","With C-lim random","With C-lim optimal"])
%fontsize(40,"points");
%xlabel("log(total metabolite concentration)");
%ylabel("log(flux efficiency)");
%title("Under C-limitation");