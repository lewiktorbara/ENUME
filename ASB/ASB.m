clear all
close all

[x10, y10] = seq(10);
[x20, y20] = seq(20);
[x30, y30] = seq(30);

F10 = lsmatrix(x10, 5);
F20 = lsmatrix(x20, 5);
F30 = lsmatrix(x30, 5);

p10 = paramiter(F10, y10);
p20 = paramiter(F20, y20);
p30 = paramiter(F30, y30);

f10 = sum(F10'.*p10);
f20 = sum(F20'.*p20);
f30 = sum(F30'.*p30);

figure('name', 'sin10')
%set(gca, 'YScale', 'log')
plot(x10,y10, 'o')
hold on
plot(x10,f10)

figure('name', 'sin20')
%set(gca, 'YScale', 'log')
plot(x20,y20, 'o')
hold on
plot(x20,f20)

figure('name', 'sin30')
%set(gca, 'YScale', 'log')
plot(x30,y30, 'o')
hold on
plot(x30,f30)

%---------------------------------------------------
i=1;
sigma = logspace(-5,-1,50);
for c=sigma
    rmsc(i,1)=1000;
    mrsc(i,1)=1000;
    for n=5:55
        for k=3:n-1
            [x,y] = seq(n);
            F = lsmatrix(x,k);
            p = paramiter(F,y);
            f = (sum(F'.*p))';
            rms(n,k) = (norm((f-y),2)/norm(y,2));
            mrs(n,k) = (norm((f-y),inf)/norm(y,inf));
            
            yc=y+randn(n,1)*c;
            pc=paramiter(F,yc);
            fc=(sum(F'.*pc))';
            
            a=norm((fc-yc),2)/norm(yc,2);
            if(a<rmsc(i,1))
                rmsc(i,1) = a;
                rmsc(i,2) = n;
                rmsc(i,3) = k;
            end
            
            a=norm((fc-yc),inf)/norm(yc,inf);
            if(a<mrsc(i,1))
                mrsc(i,1) = a;
                mrsc(i,2) = n;
                mrsc(i,3) = k;
            end
        end
    end
    i=i+1;
end
    
figure('name','rms')
surf(log10(rms));
xlabel('K');
ylabel('N');
zlabel('RMS');

figure('name','mrs')
surf(log10(mrs));
xlabel('K');
ylabel('N');
zlabel('MAX');



a = polyfit(sigma,rmsc(:,1)',1);
b = logspace(-5,-1);
c = (polyval(a,b));

figure
set(gca, 'XScale', 'log')
loglog(sigma,rmsc(:,1),'o');
hold on
loglog(b,c);


a = polyfit(sigma,mrsc(:,1)',1);
b = logspace(-5,-1);
c = (polyval(a,b));

figure
set(gca, 'XScale', 'log')
loglog(sigma,mrsc(:,1),'o');
hold on
loglog(b,c);



%-----------------------------------------------------

function [x,y] = seq(N)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

x = linspace(-1,1,N)';
y = sin(x.*pi).*exp(x- 1/3);

end

function [p] = paramiter(F, y)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

p = (transpose(F)*F)\transpose(F) * y;

end

function [F] = lsmatrix(x,K)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
F(:,1) = ones(size(x,1),1);
F(:,2) = x;
for k=3:K
    F(:,k) = 2*x.*F(:,k-1) - F(:,k-2);
end
end