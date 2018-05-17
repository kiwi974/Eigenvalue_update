% T = (diag(1:1:5) + diag(ones(4,1),1) + diag(ones(4,1),-1))
% T = (diag([0.5 0.31 0.57 0.99]) + diag([0.11 0.11 0.11],1) + diag([0.11 0.11 0.11],-1))
% no deflation necessary


function [Q,L] = divide(T)
    
    theta = pi/4;
    
	if (length(T) == 1)
		L = T;
		Q = 1;
	else 
		
		%%%%%%%%%%%%%%%%%% Constants to search zeros %%%%%%%%%%%%%%%%%%%
		itMax = 100;
		epsilon = 10^(-4);
	
		%%%%%%%%%%%%%%%%%%%%%% Partitioning of T %%%%%%%%%%%%%%%%%%%%%%%
		n = length(T);
		if (mod(n,2) == 0)
			m = n/2;
		else 
			m = (n-1)/2;
		end
		T1 = T(1:m,1:m);
		bm = T(m,m+1);
		T1(m,m) = T1(m,m) + bm;
		T2 = T(m+1:n,m+1:n);
		T2(1,1) = T2(1,1) + bm;
		u = zeros(n,1);
		u(m) = -1;
		u(m+1) = 1;
		rho = -bm;
        fill = zeros(m,n-m);
		Tconst = [T1 fill ; fill' T2] + rho*u*u';
        
		
		%%%%%%%%%%%% Call of the algorithm on T1 and T2 %%%%%%%%%%%%%%%%
		[Q1,L1] = divide(T1);  % Q1 is mxm
		[Q2,L2] = divide(T2);  % Q2 is (m+(1))x(m(+1))
		
		%%%%%%%%%%%%%%%%%%%%% Shaping of D and v %%%%%%%%%%%%%%%%%%%%%%%
		D = zeros(n,n);
		D(1:m,1:m) = diag(L1);
		D(m+1:n,m+1:n) = diag(L2);
		v = zeros(n,1);
		v(1:m) = (Q1')*u(1:m);
		v(m+1:n) = (Q2')*u(m+1:n);
		
		
		
		%%%%%%%%%%%%%%%%%%%%%%% Deflation step %%%%%%%%%%%%%%%%%%%%%%%%%
		
        [~,~,eigenvalues,eigenvectors] = deflation(D,length(D),v,rho,theta,[],[],1);
        leneig = length(eigenvalues);
        
        disp([" "])
        disp(["----------------------------------------------------------"])
        disp(["----------------------------------------------------------"])
        
        disp(["After the deflation, we got the eigenvalues : "]);
        eigenvalues
		
		
		%%%%%%%%%%%%%% Find the eigenvalues of D+rho*v*v' %%%%%%%%%%%%%%
		%%%%%%%% They are the solutions of the secular equation %%%%%%%%
		
		% Anyway the sign of rho, we always search a solution between 
		% d(i) and d(i+1) with i = 1 : n-1. Then, we have to find an 
		% additionnal zero according to the sign of rho.  
		
		lambda = [];
		
		d = diag(D); % column vector
        range = sort(diag(D));
                
        
        %figure();
        %plotsecular(@(x)secular(x,v,d,rho),-20,20,10000,d');
		
        for i = 1 : (n-1)
            alreadyFound = foundEig(leneig,eigenvalues,range(i),range(i+1),"in",rho);
            if (~alreadyFound)
                l = dichotomous(@(x)secular(x,v,d,rho),range(i),range(i+1),itMax,epsilon);
                lambda = [lambda l];
            end
        end
        
        
        disp(["For the matrix : "])
        Tconst
        disp(["Eigenvalues are : "])
        d
        disp(["And between them we found the zeros"])
        lambda
        
	    disp(["Now we have to find the last zero according to the sign of rho."])
		% Computation of the additionnal eigenvalue : we use a 
		% dichotomous search
        found = foundEig(leneig,eigenvalues,0,range(n),"out",rho);
        if ((rho > 0) && (~found))  % search between d(n) and Inf
			cs = changeSign(@(x)secular(x,v,d,rho),range(n),"i");
			lastZero = dichotomous(@(x)secular(x,v,d,rho),range(n),cs,itMax,epsilon);
			lambda = [lambda lastZero];
		else % search between -Inf and d(1) 
            if (~foundEig(leneig,eigenvalues,range(1),0,"out",rho)) 
                cs = changeSign(@(x)secular(x,v,d,rho),range(1),"d")
                lastZero = dichotomous(@(x)secular(x,v,d,rho),cs,range(1),itMax,epsilon);
                lambda = [lambda lastZero];
            end
        end
        
        lambdaQp = lambda;
        lambda = [lambda eigenvalues];
        
        disp(["All the zeros found are : "])
        lambda
        L = lambda;
        
		%%%%%%%%%%%%% Find the eigenvectors of D+rho*v*v' %%%%%%%%%%%%%%
		Qp = zeros(n,n);
        for j = 1:length(lambdaQp)
			M = (lambdaQp(j)*eye(n) - D)\v;
			Qp(:,j) = M ./ norm(M,2);
        end
        
        Qp
        
        % We complete Qp with the eigenvectors found during the deflation
        if (length(eigenvectors) >= 1)
            s = size(eigenvectors);
            for j = 1:s(1)
                Qp(:,j+length(lambdaQp)) = eigenvectors(j,:)';
            end
            Qp
        end
	
		%disp(["Computation of Qp is done and we found :"])
        %Qp 
        
		%%%%%%%%%%% Form the matrix of the eigenvectors of T %%%%%%%%%%%
		Q = zeros(n,n);
		Q(1:m,1:m) = Q1;
		Q(m+1:n,m+1:n) = Q2;
		Q = Q*Qp';
		
        disp(["Compputation of Q is done and we found :"])
        
        
        disp(["----------------------------------------------------------"])
        disp(["----------------------------------------------------------"])
        disp([" "])
		
    end	
 end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find a point where the function changes its sign, according to sense. 
% Param : -f : the function we study.
%		  -d : the departure point. 
%		  -sense : the sense of the research 
%			-> "i" : search to inverse the sign of f between d and Inf
%			-> "d" : search to inverse the sign of f between -Inf and d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = changeSign(f,di,sense)
	d = di;
	step = 100;
    if (sense == "d")
		fd = f(di-1/10000);
    else 
        fd = f(di+1/10000);
    end
	converged = false;
	while (~converged)
        if (sense == "d")
			d = d - step;
        else 
			d = d + step;
        end
        %disp([f(old)*f(d) < 0])
		converged = fd*f(d) < 0;
	end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the secular equation thanks to a dichotomous search. 
% Param : -f : the function which we search a zero.
%		  -a : the lower bound of the research interval.
%		  -b : the upper bound of the research interval.
%		  -itMax : maximum number of iterations to do. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = dichotomous(f,a,b,itMax,epsilon) 
	
    a = a + 0.00000001;
    b = b - 0.00000001;
	c = (a+b)/2;
	it = 1;

	while ((abs(f(c)) >= epsilon) && (it < itMax))
		if (f(c)*f(b) <= 0)
			a = c;
        else 
			b = c;
		end 
		c = (a+b)/2;
		it = it + 1;
	end 
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the secular function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = secular(x,v,d,rho)
	sum = 0;
	for k = 1:length(d)
		sum = sum + (v(k)^2)/(x-d(k));
	end
	y = 1 - rho*sum;
end



function [] = plotsecular(f,a,b,n,d)
    h = (b-a)/n;
    x = a:h:b;
    y = zeros(1,n+1);
    lengthd= length(d);
    for i = 1:(n+1)
        y(i)=f(x(i));
    end
    hold on
    title(strcat('Eigenvalue = ',num2str(min(d)), ' , ', num2str(max(d))))
    plot(x,y,'-b')
    plot([(a+b)/2 (a+b)/2],[1 1],'-r')
    %plot([d(lengthd) d(lengthd)],[-Inf Inf],'-m')
    hold off

end







function [f,v] = found(x,n,di,di1,place,rho)
    if (place == "in")
        f = false;
        i = 1;
        while ((i <= n) && (~f))
            f = ((di <= x) && (x <= di1));
            i = i + 1;
        end
    else
        if (rho > 0)
            f = (x > di1);
            if f 
                v = x;
            end
        else
            if (rho < 0)
                f = (x < di);
            end
        end
    end       
end




function alreadyFound = foundEig(leneig,eigenvalues,ri,ri1,place,rho)
    alreadyFound = false;
    k = 1;
    if (length(eigenvalues) >= 1)
        while ((k <= leneig) && (~alreadyFound))
            alreadyFound = found(eigenvalues(k),leneig,ri,ri1,place,rho);
            k = k + 1;
        end
    end
end