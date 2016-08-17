function ncgExample

%Graph function
[x,y] = meshgrid(-2:.2:2, -2:.2:2);                                
z = x.*exp(-x.^2 - y.^2);                                        
surf(x,y,z)


%For line search
a1 = 0;
a3 = 0.5;

%Initial point
x = [-0.6; -0.5];

for iter = 1:2
    %Calculate normalized gradient
    dz    = grad(x);
    z     = dz/norm(dz,2);
    
    %Graph 1d line along gradient
    entry = 1;
    for k = 0:1:6
        g(entry) = funct(x - k*z);
        entry    = entry+1;
    end
    figure
    plot(0:0.1:6,g)
    
    %Get alpha3 that gives a lower value than alpha1
    h1 = funct(x - a1*z);
    while funct(x - a3*z) > h1;
        a3 = a3/2;
    end

    %Check value given by alpha2 = alpha3/2
    h3 = funct(x - a3*z);
    a2 = a3/2;
    h2 = funct(x - a2*z);
    if funct(x - a2*z) > h3
        optA = a3;
    else
        optA = a2 + 1/2*( (a2-a1)^2*(h2-h3)-(a2-a3)^2*(h2-h1) )/...
                        ( (a2-a1)*(h2-h3)  -(a2-a3)*(h2-h1) );
    end
    optA
    
    %Calculate new point
    x   = x - optA*z
    val = funct(x)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = funct(x)

z = x(1,1).*exp(-x(1,1).^2 - x(2,1).^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dz = grad(x)

n = exp(-x(1,1)^2 -x(2,1)^2);

dz = n*[(1-2*x(1,1)); (-2*x(1,1)*x(2,1))];
