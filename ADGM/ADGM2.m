function [X] = ADGM2(X0, indH1, valH1, indH2, valH2, indH3, valH3, rho, MAX_ITER, verb, eta, iter1, iter2)
% input:
%  X0 / N2 x N1 matrix / initial assignment matrix
%  indHi / Nt x i matrix / indices in sparse tensor of order i
%  valHi / Nt x 1 / values in sprase tensor of order i
% output:
%  X2 / N2 x N1 matrix / output assignment matrix
    n_rhos_max = 20;
    rhos = zeros(1,n_rhos_max);
    rhos(1) = rho;
    for i=2:n_rhos_max
    	rhos(i) = rhos(i-1)*eta;
    end
    [x, ~, ~] = mexADGM_3rdORDER(X0, int32(indH1), -valH1, int32(indH2), -valH2, int32(indH3), -valH3, rhos, MAX_ITER, verb, false, iter1, iter2, 2);
    [N2,N1] = size(X0);
    X = reshape(x, N2, N1);
end
