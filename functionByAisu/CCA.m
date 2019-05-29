
function [ccaEigvector1, ccaEigvector2] = CCA(data1, data2)

 

% Input:

% data1 ¡ª¡ª view1

% data2 ¡ª¡ª view2

% both row : a sample

% column : a feature

% Output:

% ccaEigvector1 : the projection of view1

% ccaEigvector2 : the projection of view2

% both are not unit(length) one, it makes the conical

% correlation variable has unit variance

% Reference£º

% Appearance models based on kernel canonical correlation analysis

% Pattern Recognition, 2003

% Comments:

% using SVD instead of using the eigen decomposition

 

dataLen1 = size(data1, 2);

dataLen2 = size(data2, 2);

 

% Construct the scatter of each view and the scatter between them

data = [data1 data2];

covariance = cov(data);

% Sxx = covariance(1 : dataLen1, 1 : dataLen1) + eye(dataLen1) * 10^(-7);

Sxx = covariance(1 : dataLen1, 1 : dataLen1);

% Syy = covariance(dataLen1 + 1 : size(covariance, 2), dataLen1 + 1 : size(covariance, 2)) ...

% + eye(dataLen2) * 10^(-7);

Syy = covariance(dataLen1 + 1 : size(covariance, 2), dataLen1 + 1 : size(covariance, 2));

Sxy = covariance(1 : dataLen1, dataLen1 + 1 : size(covariance, 2));

% Syx = Sxy';

 

 

 

% using SVD to compute the projection

Hx = (Sxx)^(-1/2);

Hy = (Syy)^(-1/2);

 

H = Hx * Sxy * Hy;

[U, D, V] = svd(H, 'econ');

ccaEigvector1 = Hx * U;

ccaEigvector2 = Hy * V;

% make the canonical correlation variable has unit variance

ccaEigvector1 = ccaEigvector1 * diag(diag((eye(size(ccaEigvector1, 2)) ./ sqrt(ccaEigvector1' * Sxx * ccaEigvector1))));

ccaEigvector2 = ccaEigvector2 * diag(diag((eye(size(ccaEigvector2, 2)) ./ sqrt(ccaEigvector2' * Syy * ccaEigvector2))));

 

end