function ncg
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
measMu_a(rRange,cRange) = 3.0;
rRange = 5:7;
cRange = 6:8;
measMu_a(rRange,cRange) = 6.0;

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


%Graph contour plot
x = 0:0.4:10;
y = 0:0.4:10;
figure
realM = rot90(reshape(measMu_a,numPs,numPs));
[c,h] = contourf(x,y,realM);
colormap bone;
cbar_handle = colorbar('location','eastoutside');
set(get(cbar_handle,'ylabel'),'string','Absorption coefficient times c (1/ns)','FontSize',12)
daspect([1 1 1])
set(h,'EdgeColor','none') 
% xlabel('x coordinate (cm)','FontSize',16)
% ylabel('y coordinate (cm)','FontSize',16)
% title('Original phantom','FontSize',16)


%NONLINEAR CONJUGATE GRADIENT USING LINE SEARCH
%WITH PARABOLIC INTERPOLATION AND POLAK-RIVIERE
[X fX i] = minimize(optMu_a, 'getGradient', -80, Extra);


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
title('Reconstructed image from phanton using CGD (80 function evaluations)','FontSize',16)


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
title('Reconstructed image from phanton using CGD (80 function evaluations)','FontSize',16)
