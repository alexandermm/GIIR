function [x,calcRes] = bi_cgstab(d, e,f,g, h,  b, res,maxiter)
%Description:       
%Bi-STAB (BiConjugate Gradient Stabilized method with preconditionng) 
%algorithm from the paper:
%BI-CGSTAB: A fast and smoothly converging variant of 
%BI-CG for the solution of nonsymmetric linear systems
%by H.A. Van Der Vorst
%Using x0 = b./d and conditioning matrix K = diag(A) = d
%Made to handle pentadiagonal systems using vectorization


N   = length(b);                    %System size
x   = b./f;                         %Initial guess
r   = b - timesPent(d,e,f,g,h,x,N); %Get initial error r
r0  = r;

%Initialize variables
rho = 1; alpha = 1; omega = 1;
v = zeros(N,1); p = zeros(N,1);
iter = 0;
%Variables for fast multiplication

%Actual iterations
while (norm(b-timesPent(d,e,f,g,h,x,N),2)/N > res) && (iter < maxiter)
    lastRho = rho;
    rho     = r0.'*r;
    beta  = (rho/lastRho)*(alpha/omega);
    p     = r + beta*(p - omega*v);

    y     = p./f; %Can do this since K = diag(d)
    v     = timesPent(d,e,f,g,h,y,N);
    alpha = rho/(r0.'*v);
    s     = r - alpha*v;
    
    z     = s./f; %Can do this since K = diag(d)
    t     = timesPent(d,e,f,g,h,z,N); 
    Kinv  = 1./f; %Can do this since K = diag(d)
    omega = ((Kinv.*t).'*(Kinv.*s))/((Kinv.*t).'*(Kinv.*t));  
    x     = x + alpha*y + omega*z;
    
    r = s - omega*t;
    iter = iter + 1;
end
%Calculate final residual
calcRes = norm(b-timesPent(d,e,f,g,h,x,N),2)/N;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = timesPent(d, e,f,g, h,  x,N)
%Description:       
%Function for fast multiplication of pentadiagonal matrices and vectors

sqN      = sqrt(N);
numZeros = zeros(sqN,1); 
k1  = 1:N-1;
k2  = 1:N-sqN;

diags1 = [0;e.*x(k1)] + [g.*x(k1+1);0];
diags2 = [numZeros; d.*x(k2)] + [h.*x(k2+sqN); numZeros];
b = f.*x + diags1 + diags2;
