function [s,w,p,e,n] = fillA(numPs,c,D,deltat,optMu_a)

numPs2 = numPs*numPs;

s = zeros(numPs2-numPs,1);
w = zeros(numPs2-1    ,1);
p = zeros(numPs2      ,1);
e = zeros(numPs2-1    ,1);
n = zeros(numPs2-numPs,1);
%Top and bottom side points (including corner points)
for j = 1:numPs 
    p(j) = 1;
    p(j+numPs2-numPs) = 1;
end
%Right and left side points
for j = 2:numPs-1 
    p(1    +(j-1)*numPs) = 1;
    p(numPs+(j-1)*numPs) = 1;
end
%Fill rest of A (interior points)
for j = 2:numPs-1
    ind = (2:numPs-1)' + (numPs-j)*numPs;
    s(ind-numPs) = -c*D;
    w(ind-1)     = -c*D;
    p(ind  )     =  4*c*D + optMu_a(ind)*deltat + 1;
    e(ind  )     = -c*D;
    n(ind  )     = -c*D;
end