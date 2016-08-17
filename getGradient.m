function [S, dS] = getGradient(optMu_a, Extra)

numPs   = Extra.numPs;
numPs2  = Extra.numPs2;
delta   = Extra.delta;
deltat  = Extra.deltat;
D       = Extra.D;
tmax    = Extra.tmax;
impulse = Extra.impulse;
numSs   = Extra.numSs;
Y       = Extra.Y;

%Get total gradient
S  = 0;
dS = zeros(numPs2,1);
%Fill A
tic
c  = deltat/delta/delta;
[s,w,p,e,n] = fillA(numPs,c,D,deltat,optMu_a);
A = diag(s,-numPs)+diag(w,-1)+diag(p)+diag(e,1)+diag(n,numPs);
invtA = inv(A)';

for source = 1:numSs 
    [Sen, dSdM] = AD(numPs,deltat,tmax,impulse,source,numSs,Y(:,source),invtA);
    S  = S  + Sen; 
    dS = dS + dSdM;
end
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,dSdM] = AD(numPs,deltat,tmax,impulse,source,numSs,Yn,invtA)


%GET Un

%Constants
numPs2  = numPs*numPs;

%Solution matrix
U = zeros(numPs2,1);

%What receivers are being used?
receivers = 1:numSs;
receivers(source) = [];

%Get source index (in natural ordering) and apply initial condition
sInd = getSourceInd(source);
U(sInd) = impulse;

%Get receiver indices
rInd = getReceiverInd(receivers);

%Run ADI scheme
tnum = 0;
Un = zeros(numPs2,tmax);
while tnum < tmax
    U = invtA'*U;
            
    %Store info
    Un(:,tnum+1) = U;
    
    %Increase timestep
    tnum = tnum + 1;
end

%GET S and delSdelU
numRs = numSs-1;
change   = zeros(numPs2,tmax);
delSdelU = zeros(numPs2,tmax);
for n = 1:tmax
    change(rInd,n)   = (Un(rInd,n) - Yn((1:numRs)'+(n-1)*numRs)).^2;
    delSdelU(rInd,n) =  Un(rInd,n) - Yn((1:numRs)'+(n-1)*numRs);  
end
S = sum(sum(change))/2;

dSdU(:,tmax) = delSdelU(:,tmax); 
for n = (tmax-1):-1:1
    %Notice the transpose of A is being used
    dSdU(:,n) = invtA*dSdU(:,n+1) + delSdelU(:,n);
end


%GET intPs and totalNPs
%These vector and value are used to determine if a diagonal point is 
%in the interior of the real matrix 
intPs = zeros(numPs2,1);
for j = 2:numPs-1
    ind = (2:numPs-1)' + (numPs-j)*numPs;
    intPs(ind) = 1;
end

%GET dSdM
dSdM = zeros(numPs2,1);
for n = 1:tmax
    %Get X for given timestep
    X = diag(intPs.*(-deltat));   
    X = X.*repmat(Un(:,n),1,numPs2);
    
    %Accumulate result
    dSdM = dSdM + X'*invtA*dSdU(:,n); 
end