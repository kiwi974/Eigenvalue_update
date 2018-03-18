function [Q,L] = arbenz(T)

	if (length(T) == 1)
		L = T;
		Q = 1;
	else 
		
		%%%%%%%%%%%%%%%%%% Constants to search zeros %%%%%%%%%%%%%%%%%%%
		itMax = 100;
		epsilon = 10^(-3);
	
		%%%%%%%%%%%%%%%%%%%%%% Partitioning of T %%%%%%%%%%%%%%%%%%%%%%%
		n = length(T);
		if (mod(n,2) == 0)
			m = n/2;
		else 
			m = (n-1)/2;
		end
		T1 = T(1:m,1:m);
		bm = T(m,m+1);
		T1(m,m) = T1(m,m) + bm
		T2 = T(m+1:n,m+1:n);
		T2(1,1) = T2(1,1) + bm
		u = zeros(n,1);
		u(m) = -1;
		u(m+1) = 1;
		rho = -bm;
		% T = [T1 0 ; 0 T2] + rho*u*u'
		
		
		%%%%%%%%%%%% Call of the algorithm on T1 and T2 %%%%%%%%%%%%%%%%
		[Q1,L1] = arbenz(T1);  % Q1 is mxm
		[Q2,L2] = arbenz(T2);  % Q2 is (m+(1))x(m(+1))
		
		
		%%%%%%%%%%%%%%%%%%%%% Shaping of D and v %%%%%%%%%%%%%%%%%%%%%%%
		D = zeros(n,n);
		D(1:m,1:m) = L1;
		D(m+1:n,m+1,n) = L2;
		v = zeros(n,1);
		v(1:m) = Q1'*u(1:m);
		v(m+1:n) = Q2'*u(m+1,n);
		
		
		%%%%%%%%%%%%%%%%%%%%%%% Deflation step %%%%%%%%%%%%%%%%%%%%%%%%%
		%%%?????????????????????????????????????????????????????????????
		%%%????????????????????????? TO DO ?????????????????????????????
		%%%?????????????????????????????????????????????????????????????
		%%%?????????????????????????????????????????????????????????????
		
		
		%%%%%%%%%%%%%% Find the eigenvalues of D+rho*v*v' %%%%%%%%%%%%%%
		%%%%%%%% They are the solutions of the secular equation %%%%%%%%
		
		% Anyway the sign of rho, we always search a solution between 
		% d(i) and d(i+1) with i = 1 : n-1. Then, we have to find an 
		% additionnal zero according to the sign of rho.  
		
		L = zeros(n-1,1);  % then we will concatenant the last zero
		
		d = diag(D); % column vector
		
		for i = 1 : (n-1) 
			L(i) = dichotomous(@(x)secular(x,v,d),d(i),d(i+1),itMax,epsilon);
		end 
	
		% Computation of the additionnal eigenvalue : we use a 
		% dichotomous search
		lastZero = 0;
		if (rho > 0) % search between d(n) and Inf
			cs = changeSign(@(x)secular(x,v,d),d(n),"i");
			lastZero = dichotomous(@(x)secular(x,v,d),d(n),cs,itMax,epsilon);
		else % search between -Inf and d(1) 
			cs = changeSign(@(x)secular(x,v,d),d(1),"d");
			lastZero = dichotomous(@(x)secular(x,v,d),cs,d(1),itMax,epsilon);
		end
		
		%%%%%%%%%%%%% Find the eigenvectors of D+rho*v*v' %%%%%%%%%%%%%%
		Qp = zeros(n,n);
		
		
		%%%%%%%%%%% Form the matrix of the eigenvectors of T %%%%%%%%%%%
		Q = zeros(n,n);
		Q(1:m,1:m) = Q1;
		Q(m+1:n,m+1:n) = Q2;
		Q = Q*Qp';
		
		
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

function x = changeSign(f,d,sense):
	step = 50;
	converged = false;
	while (~converged)
		old = d;
		if (sense == "d")
			d = d - step;
		else 
			d = d + step;
		end
		converged = f(old)*f(d) < 0;
	end
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the secular function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = secular(x,v,d)
	sum = 0;
	for k = 1:length(d)
		sum = sum + v(k)^2/(x-d(k));
	end
	y = 1 - rho*sum;
end
