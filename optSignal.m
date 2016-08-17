function signal = optSignal(numPs,delta,deltat,tmax,impulse,source,numSources,D,Mu_a)
%Name:           Alejandro P. Martinez
%Language used:  MATLAB language
%Description:       
%Provides 15 vectors giving the return signal by solving the 
%time-dependent diffusion equation for photon propagation in 
%a turbid medium using alternating direction implicit (ADI) 
%finite differences with zero Diritchlet boundary conditions

%Bi-cgstab variables
res     = 1e-12;
maxiter = 100;

%Obtaining numPs2
numPs2 = numPs*numPs;
%Constants
c      = deltat/delta/delta;

%Solution matrix
U = zeros(numPs2,1);

%What receivers are being used?
receivers = 1:numSources;
receivers(source) = [];

%Get source index (in natural ordering) and apply initial condition
sInd = getSourceInd(source);
U(sInd) = impulse;

%Display initial conditions
%realU = flipud(rot90(reshape(Uzero,numPs,numPs)));
%figure
%surf(realU)

%Get receiver indices
rInd = getReceiverInd(receivers);

%Fill A
[s,w,p,e,n] = fillA(numPs,c,D,deltat,Mu_a);

%Run ADI scheme
tnum = 0;
numRs = numSources-1;
signal  = zeros(tmax*numRs,1);
while tnum < tmax
    U = bi_cgstab(s, w,p,e, n,  U, res,maxiter);
    
    %Display data
%     if mod(tnum,2) == 0
%         figure
%         realU = rot90(reshape(U,numPs,numPs));
%         surf(flipud(realU))
%     end
    
    %Store info
    signal((1:numRs)'+tnum*numRs) = U(rInd);
    
    %Increase timestep
    tnum = tnum + 1;
end