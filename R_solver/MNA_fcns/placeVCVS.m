function [B,A] = placeVCVS(B,A,VCVS_nodes,gain)
%VCVS nodes:
% 1: positive voltage into VCVS, positive node of Rin (alpha)
alpha = VCVS_nodes(1);
% 2: negative voltage into VCVS, negative node of Rin (beta)
beta = VCVS_nodes(2);
% 3: positive node of VCVS, negative node of Rout (gamma)
gamma = VCVS_nodes(3);
% 4: negative node of VCVS (delta)
delta = VCVS_nodes(4);

B(alpha) = plus(B(alpha),-gain);
B(beta) = plus(B(beta),gain);
B(gamma) = plus(B(gamma),1);
B(delta) = plus(B(delta),-1);

A(gamma) = plus(A(gamma),1);
A(delta) = plus(A(delta),-1);

end