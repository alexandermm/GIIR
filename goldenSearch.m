function goldenSearch
%Test nonlinear conjugate gradient


clear all;
format long;

%GEOMETRIC PROPERTIES
delta   = 0.4;  %cm
numPs   = 26;
numPs2  = numPs*numPs;
deltat  = 1;    %ns
%Number of time steps for simulating signal
tmax    = 10;
%INITIAL CONDITIONS
impulse = 100;  %W/cm^2

%Number of sources
numSs = 16;

%MATERIAL PROPERTIES
%Constant difussion
D       = 2.0; %cm^2/ns
%Background mu_a
Bmu_a   = 1.0; %1/ns
%Measured mu_a
measMu_a = repmat(Bmu_a, numPs,numPs);%readBrain(numPs);
%Mu_a: vector that is optimized
optMu_a  = repmat(Bmu_a, numPs2,1);
%Object row and column ranges
rRange = 11:15;
cRange = 18:22;
measMu_a(rRange,cRange) = 10.0;
rRange = 5:7;
cRange = 6:8;
measMu_a(rRange,cRange) = 15.0;

%Display object
%figure
%surf(flipud(measMu_a))

%Change matrix to vector
measMu_a = reshape(rot90(measMu_a,-1),numPs2,1);

%Get Y
for source = 1:numSs 
    Y(:,source) = ...
    optSignal(numPs,delta,deltat,tmax,impulse,source,numSs,D,measMu_a);
end


%Get total gradient
dS = zeros(numPs2,1);
%Fill A
c  = deltat/delta/delta;
[s,w,p,e,n] = fillA(numPs,c,D,deltat,optMu_a);
A = diag(s,-numPs)+diag(w,-1)+diag(p)+diag(e,1)+diag(n,numPs);
invtA = inv(A)';

for source = 1:numSs 
    dSdM = AD(numPs,deltat,tmax,impulse,source,numSs,Y(:,source),invtA);
    dS = dS + dSdM;
end

%Golden search
changeM = 0.1;
deltaM  = changeM/max(abs(dS));
gRatio  = 1.618;
for step = 1:10
    %Calculate new S
    %Get U
    for source = 1:numSs 
        U(:,source) = ...
        optSignal(numPs,delta,deltat,tmax,impulse,source,numSs,D,optMu_a);
    end
    S = sum(sum((U-Y).^2))/2
    
    %Calculate new point
    optMu_a = optMu_a - dS*gRatio*deltaM;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dSdM = AD(numPs,deltat,tmax,impulse,source,numSs,Yn,invtA)


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
delSdelU = zeros(numPs2,tmax);
for n = 1:tmax
    delSdelU(rInd,n) =  Un(rInd,n) - Yn((1:numRs)'+(n-1)*numRs);  
end

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