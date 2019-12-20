
function [a , b] = GradNusMod( X,Y,Z )
clc
    r=100;
    ydaug =0.001;
    iter=0;
    ifun = 0;
    Xi11 = X;
    Xi22 = Y;
    Xi33 = Z;
    %a=2;b=2;c=2;
    pointsListX = [X];
    pointsListY = [Y];
    pointsListZ = [Z];
    ListSize = 1;
    e = 1/10000;
    while (abs(X+Y+Z - 1) > 0.005) %a > e && b > e && c > e%(Gi(Xi11,Xi22,Xi33) ~= 1)
        gi = Gi(X,Y,Z);
        Xi1 = Xi11 - ydaug * GradX(Xi11,Xi22,Xi33,r);
        Xi2 = Xi22 - ydaug * GradY(Xi11,Xi22,Xi33,r);
        Xi3 = Xi33 - ydaug * GradZ(Xi11,Xi22,Xi33,r);
        %a=abs(GradX(Xi1,Xi22,Xi33,r)-GradX(Xi11,Xi22,Xi33,r));
        %b=abs(GradZ(Xi1,Xi22,Xi33,r)-GradZ(Xi11,Xi22,Xi33,r));
        %c=abs(GradY(Xi1,Xi22,Xi33,r)-GradY(Xi11,Xi22,Xi33,r));
        %disp('------------Iteracija----')
        %disp(GradX(Xi11,Xi22,Xi33,r))
        %disp(GradY(Xi11,Xi22,Xi33,r))
        %disp(GradZ(Xi11,Xi22,Xi33,r))
        %disp('----------------')
        ifun = ifun + 3;
        Xi11 = Xi1
        Xi22 = Xi2
        Xi33 = Xi3
        pointsListX = [pointsListX, Xi11];
        pointsListY = [pointsListY, Xi22];
        pointsListZ = [pointsListZ, Xi33];
        r = r * 0.95 ;
        ListSize = ListSize + 1;
 %disp('Iteracija')
 %disp(i)
        %disp(Xi1)
%disp(Xi2)
%disp(Xi3)
 disp('---------------------------------')
        iter=iter+1;
    end
    
    Xt = 0:0.01:1;
    Yt = 0:0.01:1;
    Zt = 0:0.01:1;
    size = 100;
    Xte = [];
    Yte = [];
    Zte = [];
    eSize = 1000000;
    for i = 1:size
        for j = 1:size
            Xte = [Xte, Xt(i)];
            Yte = [Yte, Yt(j)];
            Zte = [Zte, 1-Xt(i)-Yt(j)];
        end
    end
    
    plot3(pointsListX,pointsListY,pointsListZ,'o',Xte,Yte,Zte,'.');
    a = Xi11;
    b = Xi22;
    fprintf('iteracij� skai�ius: %d\n',iter);
    fprintf('funkcij� skai�iavim� skai�ius: %d\n',ifun);
    fprintf('Gauti sprendiniai: %12.14f %12.14f %12.14f\n',Xi11,Xi22,Xi33);
    fprintf('Funkcijos minimumo �vertis: %12.14f',Function(Xi11,Xi22,Xi33));
end

function x = GradX(X,Y,Z,r)
    x = -1/8 * Y * Z + 2*r * (X+Y+Z-1)
end

function y = GradY(X,Y,Z,r)
    y = -1/8 * X * Z + 2*r * (X+Y+Z-1)
end

function z = GradZ(X,Y,Z,r)
    z = -1/8 * Y * X + 2*r * (X+Y+Z-1)
end

function z = Function(X,Y,Z)
    z = -1/8.* X * Y * Z;
end


function z = BaudFunction(X,Y,Z,r)
    z = -1/8.* X * Y * Z + r * ((X+Y+Z-1)^2);
end

function g = Gi(X,Y,Z)
  if abs(X+Y+Z - 1)<0.005
      g = 1;
  elseif X+Y+Z ~= 1
      g = 0; 
  end
end