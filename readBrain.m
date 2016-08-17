function img = readBrain(numPs)


%Get image
img = imread('brain','jpg');
%Graph image
% image(img)
% colormap bone
% pbaspect([1 1 1])

%Change image size
size = numPs - 2;
img = imresize(img,[size size],'bicubic');

%Graph modified image
% figure
% image(img)
% colormap bone
% pbaspect([1 1 1])


%Change values to match needed cMu 
img = double(img);

%Fix corners
%Not exactly equal to 0; might be due to bicubic interpolation 
corners = [1:2 size-2:size]; 
img(corners,corners) = 0.0;

maxColor = max(max(img));
minColor = min(min(img));

maxMua = 10.0;
minMua =  1.0;

%change the values
img = img.*(maxMua-minMua)/(maxColor-minColor) + minMua;

%Put edges
back = repmat(minMua,numPs,numPs);
back(2:numPs-1,2:numPs-1) = img;
img = back;

%Show surface
% figure
% surf(0:numPs-1,0:numPs-1,img)
% colormap bone
% shading interp



