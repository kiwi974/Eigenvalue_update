function [D,v,v_prime,eigenvalues,eigenvectors, n_deflated, G] = deflation2(d,n,v)
% Rotation matrix
DD = d;  
D = d;
d = diag(d);
G = eye(n);
tol = 1e-6;
for j = 1: n-1
    if ( abs(d(j) - d(j+1)) < tol)
        [c ,s] = rotate(v(j),v(j+1));
        G_temp = eye(n);
        G_temp(j,j) = c; G_temp(j,j+1) = s; G_temp(j+1,j) = -s; G_temp(j+1,j+1) = c;
        G =G_temp*G;
    end
end

% Deflation - simplify D and v
n_deflated = 0;
index_deflation = zeros(0);
v = G'*v; % rotate v
v_prime = zeros(n_deflated,1)'; % v without zeros
kk = 1;
for i = 1:n
    if abs( v(i)) < tol
         index_deflation(length(index_deflation)+1) = i;
    else
        n_deflated = n_deflated +1;
        v_prime(kk) = v(i);
        kk = kk+1;
    end  
end
ii = 1;
nn = n;
i = 1;
% Remove entries of D corresponding to deflated entries 
while (ii <= nn)
 if abs( v(i)) < tol
     D(ii, :) = [];
     D(:, ii) = [];
     nn = nn-1;
 else
    ii = ii+1;
 end
i = i+1;
end

eigenvectors = zeros(n,n);
eigenvalues = zeros(n,n); 
% Note that these will contain the final eigenvalues and eigenvectors (dimension is already n x n)
% Compute eigenpairs for deflated cases
ii = 1;
for i=1:n
     if abs( v(i)) < tol
        % Eigenpair for deflated case
        eigenvalues(i,i) = DD(index_deflation(ii),index_deflation(ii));
        ii = ii+1;
        eigenvectors(:,i) = e(i,n);
     end
end
end






function [c,s] = rotate(x1,x2)
% Compute Given`s rotation parameters
x = [abs(x1) abs(x2)];
c = x1/norm(x);
s = -x2/norm(x);
end


function ei = e(i,n)
ei = zeros(n, 1);
ei(i) = 1;
end