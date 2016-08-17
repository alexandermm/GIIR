function sd
%Test nonlinear gradient descent with AD


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

%Put all constants in Extra structure
Extra.numPs   = numPs;
Extra.numPs2  = numPs2;
Extra.delta   = delta;
Extra.deltat  = deltat;
Extra.D       = D;
Extra.tmax    = tmax;
Extra.impulse = impulse;
Extra.numSs   = numSs;
Extra.Y       = Y;












%NONLINEAR CONJUGATE GRADIENT USING SIMPLE LINE SEARCH
%WITH PARABOLIC INTERPOLATION AND STEEPEST DESCENT
a1 = 0;
a3 = 5;
numEvals = 0;
place = 1;
while numEvals < 80
    [h1, ds] = getGradient(optMu_a, Extra);
    numEvals = numEvals+1
    
    %Record information
    sVal(place) = h1;
    eVal(place) = numEvals;
    place = place+1;
    
    
    %Calculate normalized direction 
    z = -ds/norm(ds,2);
    
    %Graph 1d line along gradient
%     entry = 1;
%     for k = 0:1:6
%         g(entry) = getGradient(optMu_a + k*z, Extra);
%         entry    = entry+1;
%     end
%     figure
%     plot(0:1:6,g)
    
    %Get alpha3 that gives a lower value than alpha1
    numDivs = 0;
    numEvals = numEvals+1;
    while (getGradient(optMu_a + a3*z, Extra) > h1) && (numDivs < 10)
        a3 = a3/2;
        numDivs = numDivs+1;
        numEvals = numEvals+1;
    end

    %Check value given by alpha2 = alpha3/2
    h3 = getGradient(optMu_a + a3*z, Extra);
    numEvals = numEvals+1;
    a2 = a3/2;
    h2 = getGradient(optMu_a + a2*z, Extra);
    numEvals = numEvals+1;
    if h2 > h3
        optA = a3;
    else
        optA = a2 + 1/2*( (a2-a1)^2*(h2-h3)-(a2-a3)^2*(h2-h1) )/...
                        ( (a2-a1)*(h2-h3)  -(a2-a3)*(h2-h1) );
    end
    optA
    %Calculate new point
    optMu_a = optMu_a + optA*z;
    
    %Graph new mua distribution
%     if mod(iter,1) == 0
%         realM = rot90(reshape(optMu_a,numPs,numPs));
%         figure
%         surf(flipud(realM))
%     end
end


figure
plot(eVal,sVal)
xlabel('Function evaluation','FontSize',16)
ylabel('Value of sensitivity (1/ns^2)','FontSize',16)
title('Value of sensitivity vs. number of function evaluations (SD)','FontSize',16)

x = 0:0.4:10;
y = 0:0.4:10;
X = optMu_a;
%Graph surface plot
realM = rot90(reshape(X,numPs,numPs));
figure
surf(x,y,realM)
daspect([1 1 1])
colormap bone
shading  interp
xlabel('x coordinate (cm)','FontSize',16)
ylabel('y coordinate (cm)','FontSize',16)
zlabel('Absorption coefficient times c (1/ns)','FontSize',16)
title('Reconstructed image from phanton using SD (80 function evaluations)','FontSize',16)


%Graph contour plot
figure
[c,h] = contourf(x,y,realM);
colormap bone;
cbar_handle = colorbar('location','eastoutside');
set(get(cbar_handle,'ylabel'),'string','Absorption coefficient times c (1/ns)','FontSize',12)
daspect([1 1 1])
set(h,'EdgeColor','none') 
xlabel('x coordinate (cm)','FontSize',16)
ylabel('y coordinate (cm)','FontSize',16)
title('Reconstructed image from phanton using SD (80 function evaluations)','FontSize',16)














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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