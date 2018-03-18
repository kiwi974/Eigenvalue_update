function [] = solver(num) 
	if (num == 2)
		test1()
	else 
		disp(['No solver for this']);
	end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the secular equation with adequated Newtons method.
% Param : -v : the vector such that v*v^T it ourmatrix perturbation.
%		  -d : the vector of the eigenvalues of D
%		  -i : the indice of the interval in which we want to find a zero 
%		  -l0 : an approximation of the zero of the secular function in 
%				[d(i),d(i+1)]
%		  -pho : in that first version, we assume that pho = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [l] = secularSolver(v,d,i,l0,tol)

	n = length(v);

	% Compute c1, c2, c3
	c1 = dphi1(l0,v,d,i)*(d(i)-l0)^2;
	cc1 = phi1(l0,v,d,i) - c1;

	c2 = dphi2(l0,v,d,i,n)*(d(i+1)-l0)^2;
	cc2 = phi2(l0,v,d,i,n) - c2;
	
	c3 = 1 + cc1 + cc2;
	
	% Solve h(x) = 0
	l = newton(@(x)h(x,c1,c2,c3,d,i));
end




%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define the function h that locally belongs to f
function y = h(x,c1,c2,c3,d,i)
	y = c3 + c1/(d(i)-x) + c2/(d(i+1)-x);
end

% y = phi1(x)
function y = phi1(x,v,d,i)
	y = 0
	for k = 1:i
		y = y + (v(k)^2)/(d(k)-x);
	end
end 


%y = phi1'(x)
function y = dphi1(x,v,d,i)
	y = 0
	for k = 1:i
		y = y + x*(v(k)^2)/((d(k)-x)^2);
	end
end


% y = phi2(x)
function y = phi2(x,v,d,i,n)
	y = 0
	for k = (i+1):n
		y = y + (v(k)^2)/(d(k)-x);
	end
end 


%y = phi2'(x)
function y = dphi2(x,v,d,i,n)
	y = 0
	for k = (i+1):n
		y = y + x*(v(k)^2)/((d(k)-x)^2);
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
	
	if (it >= itMax)
		flag = 2
	else 
		flag = 1
	end
end



% Test with simple functions 
function test1()
	function y = s(x)
		y = x^3 + 6*x^2 - 5;
	end
	zero = dichotomous(@s,-4,3,100,10^(-3)) 
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

