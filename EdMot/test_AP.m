%test ap : 2018-12-20
function results = test_AP(A,ty_lcc)
[i,j,s] = find(A);
S = [i,j,s];
p = median(s);
MAX_iter = 1;
nmi_tmp = zeros(MAX_iter,1);
pty_tmp = zeros(MAX_iter,1);
f1_tmp = zeros(MAX_iter,1);
Q_tmp = zeros(MAX_iter,1);
for i=1:MAX_iter
    [com,netsim,dpsim,expref]=apcluster(S,p);
     nmi_tmp(i) = NMI(ty_lcc,com);
     pty_tmp(i) = Purity(ty_lcc,com);
     f1_tmp(i) = F1Over(ty_lcc,com);
     Q_tmp(i) = modul(A,com);
end

 nmi = mean(nmi_tmp); std_nmi = std(nmi_tmp);
 pty = mean(pty_tmp); std_pty = std(pty_tmp);
 f1 = mean(f1_tmp); std_f1 = std(f1_tmp);
 Q = mean(Q_tmp); std_Q = std(Q_tmp);
 results = [nmi,std_nmi,pty,std_pty,f1,std_f1,Q,std_Q];
end
