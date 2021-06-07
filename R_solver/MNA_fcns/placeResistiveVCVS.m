function [Y,B,A,D22] = placeResistiveVCVS(Y,B,A,D22,VCVS_nodes,gain,Rin,Rout,numEtc)

%VCVS nodes:
% 1: positive voltage into VCVS, positive node of Rin (alpha)
alpha = VCVS_nodes(1);
% 2: negative voltage into VCVS, negative node of Rin (beta)
beta = VCVS_nodes(2);
% 3: positive node of VCVS, negative node of Rout (gamma)
gamma = VCVS_nodes(3);
% 4: negative node of VCVS (delta)
delta = VCVS_nodes(4);

Y(alpha,alpha) = plus(Y(alpha,alpha),1/Rin);
Y(alpha,beta) = plus(Y(alpha,beta),-1/Rin);
Y(beta,alpha) = plus(Y(beta,alpha),-1/Rin);
Y(beta,beta) = plus(Y(beta,beta),1/Rin);

B(alpha) = plus(B(alpha),-gain);
B(beta) = plus(B(beta),gain);
B(gamma) = plus(B(gamma),1);
B(delta) = plus(B(delta),-1);

A(gamma) = plus(A(gamma),1);
A(delta) = plus(A(delta),-1);

D22(numEtc,numEtc) = plus(D22(numEtc,numEtc),Rout);

end