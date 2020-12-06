I = imread('door01.jpg');
[Hs W col]=size(I);
%imshow(I)
% [x,y] = getpts;
% xp1=x(1,1);
% xp2=x(2,1);
% xp3=x(3,1);
% xp4=x(4,1);
%Points plucked out for four sides of frame
xp1 = 186;
xp2=369;
xp3=371;
xp4=192;


% yp1=y(1,1);
% yp2=y(2,1);
% yp3=y(3,1);
% yp4=y(4,1);
yp1=101;
yp2=119;
yp3=303;
yp4=332;


%For Frame
x1=1;
x2=77.5;
x3=77.5;
x4=1;

y1=1;
y2=1;
y3=62;
y4=62;

%For Door
% x1=0;
% x2=76.5;
% x3=76.5;
% x4=0;
% 
% y1=0;
% y2=0;
% y3=61;
% y4=61;

A=[
-x1  -y1  -1   0    0    0   x1*xp1   y1*xp1   xp1;
 0    0    0 -x1   -y1  -1   x1*yp1   y1*yp1   yp1;
-x2  -y2  -1   0    0    0   x2*xp2   y2*xp2   xp2;
 0    0    0 -x2   -y2  -1   x2*yp2   y2*yp2   yp2;
-x3  -y3  -1   0    0    0   x3*xp3   y3*xp3   xp3;
 0    0    0 -x3   -y3  -1   x3*yp3   y3*yp3   yp3;
-x4  -y4   -1  0    0    0   x4*xp4   y4*xp4   xp4;
 0    0    0  -x4  -y4  -1   x4*yp4   y4*yp4   yp4];

[U,S,V] = svd(A);


H=V(:,end);
H=reshape(H,3,3);


invH = inv(H);
projected_point1=[1, 1, 1]*invH;
xn=projected_point1(1,1)/projected_point1(1,3);
yn=projected_point1(1,2)/projected_point1(1,3);

projected_point2=[W, 1, 1] * invH;
xn2=projected_point2(1,1)/projected_point2(1,3);
yn2=projected_point2(1,2)/projected_point2(1,3);

projected_point3=[W, Hs, 1] * invH;
xn3=projected_point3(1,1)/projected_point3(1,3);
yn3=projected_point3(1,2)/projected_point3(1,3);

projected_point4=[1, Hs ,1] * invH;
xn4=projected_point4(1,1)/projected_point4(1,3);
yn4=projected_point4(1,2)/projected_point4(1,3);
xl=xn3-xn4;
yl=yn3-yn2;
xl=round(xl);
yl=round(yl);

x=zeros(yl,xl,3,'uint8');

%Reverse
% for i=1:800
% for j=1:600
% projected_point=[i j 1] * invH;
% xi=projected_point(1,1)/projected_point(1,3);
% yi=projected_point(1,2)/projected_point(1,3);
% for k=1:3
% x(floor(yi+128),floor(xi+222),k)=I(j,i,k);
% 
% end
% end
% end
%imshow(x)
%projected_point=[1 1 1]*invH
% xi=projected_point(1,1)/projected_point(1,3)
% yi=projected_point(1,2)/projected_point(1,3)
xmax=0;
ymax=0;
for i=1:xl %753
for j=1:yl %527
projected_point=[i-55, j-130, 1]*H;
xi=projected_point(1,1)/projected_point(1,3);
yi=projected_point(1,2)/projected_point(1,3);
xc=floor(xi);
yc=floor(yi);
xmax(xc>xmax)=xc;
ymax(yc>ymax)=yc;
xc(xc<1)=0;
yc(yc<1)=0;
xc(xc>800)=0;
yc(yc>600)=0;
for k=1:3
if(xc~=0 && yc~=0)
 x(j,i,k)=I(yc,xc,k);
else
 x(j,i,k)=0;
end

end
end
end
% imwrite(x,'x.jpg');
%imshow(x)


%imwrite(img2,'ex.jpg');

% Maintain Width
ratio=W/xl;
H1 = round(yl*ratio);

Y = zeros(H1,W,3,'uint8');
hs = (yl/H1)
ws = (xl/W)

for k=1:3
    for i=1:H1-1
      yic = (hs * i) + (0.5 * (1 - 1/ratio));
       for j=1:W-1
           xic = (ws * j) + (0.5 * (1 - 1/ratio));
      %// Any values out of acceptable range
          xic(xic < 1) = 1;
          xic(xic > W - 0.001) = W - 0.001;
          xi1 = floor(xic);
          xi2 = xi1 + 1;
          yic(yic < 1) = 1;
          yic(yic > xl - 0.001) = xl - 0.001;
          yi1 = floor(yic);
          yi2 = yi1 + 1;
      %// 4 Neighboring Pixels
          NP1 = x(yi1,xi1,k);
          NP2 = x(yi1,xi2,k);
          NP3 = x(yi2,xi1,k); 
          NP4 = x(yi2,xi2,k);
      %// 4 Pixels Weights
          PW1 = (yi2-yic)*(xi2-xic);
          PW2 = (yi2-yic)*(xic-xi1);
          PW3 = (xi2-xic)*(yic-yi1);
          PW4 = (yic-yi1)*(xic-xi1);
          Y(i,j,k) = PW1 * NP1 + PW2 * NP2 + PW3 * NP3 + PW4 * NP4;
        end
    end
end

imshow(Y)
imwrite(Y,'doorresult.jpg');


%Test
% tform = projective2d(A); 
% img2 = imwarp(I, tform); 
% imshow(img2)
% imresize
% end
% ratiox=size(x,1)/600
% ratioy=size(x,2)/800
% for h=1:3
%     s(:,:,h)=trad_bil(I,ratio1,ratio2)
% I(1,1,3)
% x(34,245,1)
