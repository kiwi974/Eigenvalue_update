%% Deflation step 

function [Dnd,vnd,eigenvalues,eigenvectors] = deflation(D,n,v,rho,theta,eigenvalues,eigenvectors,ibegin)

    Dnd = D;
    vnd = v;
    
    C = 10e-3;
    
    normT = norm(D+rho*(v*v'),2);
    
    for i = ibegin:n
        normvi = norm(v(i),1);
        if (normvi < C*normT)
            eigenvalues = [eigenvalues D(i,i)];
            eigenvectors = [eigenvectors ; e(i,n)];
        else 
            for j = (i+1):n
                if (abs(D(i,i)-D(j,j)) < C*normT)
                    G = givens(n,i,j,theta);
                    [Dnd,vnd,eigenvalues,eigenvectors] = deflation(G'*D,n,G'*v,rho,theta,eigenvalues,eigenvectors,j);
                end
            end
        end
    end   
end






function G = givens(n,i,j,theta)
    G = diag(ones(n,1));
    G(i,i) = cos(theta);
    G(j,j) = cos(theta);
    G(j,i) = sin(theta);
    G(i,j) = -sin(theta);
end


function ei = e(i,n)
    ei = zeros(1,n);
    ei(i) = 1;
end



%% Testing 

%D = diag([0.1981 1.5550 3.2470 0.1981 1.5550 3.2470]);
%v = [0.7370 -0.5910 0.3280 0.7370 -0.5910 0.3280]';
%rho = -1;
%theta = pi/4;
%eigenvalues = [];
%eigenvectors = [];
%[Dnd,vnd,eigenvalues,eigenvectors] = deflation(D,6,v,rho,theta,eigenvalues,eigenvectors,1)