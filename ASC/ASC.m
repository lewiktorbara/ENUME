clear
close all

A = [0, 1;-1.25 , -1];
h = 0.01;
y0 = [4; -2];
ye(:,1)  =  [4 -2];
u(:,1) = [4 -2];

I = eye(2);
Bu = [I-(5/36)*h*A, ((sqrt(15)/15)-(2/9))*h*A, ((sqrt(15)/30)-5/36)*h*A;
      -(5/36 + sqrt(15)/24)*h*A, I-(2/9)*h*A, ((sqrt(15)/24)-5/36)*h*A; 
      -(5/36 + sqrt(15)/30)*h*A, -(2/9 + sqrt(15)/15)*h*A, I-(5/36)*h*A];

  
opts = odeset('Reltol',10^(-13),'AbsTol',10^(-14),'Stats' ,'on');
t = 0:h:10;
[ode_t, ode_y] = ode113(@fun,[0 10],y0',opts);


  

for i=2:length(t)
    F(:,i) = Bu\[A*u(:,i-1);A*u(:,i-1) ;A*u(:,i-1)];
    u(:,i) = u(:,i-1) + h*((5/18)*F(1:2,i) + (4/9)*F(3:4,i) +(5/18)*F(5:6,i));
        
    ye(:,i) = ye(:,i-1) + A*h*ye(:,i-1);
    
end

figure('name','1')
plot(t,u(1,:),'-o');
hold on
plot(t,ye(1,:),'-o');
hold on
plot(ode_t, ode_y(:,1),'-');
hold off
legend('GAUSS-LEGENDRE','EULER','ODE113')





%------------------------------------------------------------   
hspace = [logspace(-5,0,20) 2 3 4 5];
iter = 1;    
for h=hspace
    clear u1 ode1_y t1 ode1_t F ye1;
    
    u1(:,1) = [4 -2];
    ye1(:,1) = [4 -2];
    Bu = [I-(5/36)*h*A, ((sqrt(15)/15)-(2/9))*h*A, ((sqrt(15)/30)-5/36)*h*A;
      -(5/36 + sqrt(15)/24)*h*A, I-(2/9)*h*A, ((sqrt(15)/24)-5/36)*h*A; 
      -(5/36 + sqrt(15)/30)*h*A, -(2/9 + sqrt(15)/15)*h*A, I-(5/36)*h*A];


    t1 = 0:h:10;
    [ode1_t, ode1_y] = ode113(@fun,t1,y0',opts);


 
    for i=2:length(t1)
        F(:,i) = Bu\[A*u1(:,i-1);A*u1(:,i-1) ;A*u1(:,i-1)];
        u1(:,i) = u1(:,i-1) + h*((5/18)*F(1:2,i) + (4/9)*F(3:4,i) +(5/18)*F(5:6,i));
        
        ye1(:,i) = ye1(:,i-1) + A*h*ye1(:,i-1);
        
    end
      
    rmse(iter) = norm(ye1(1,:)'-ode1_y(:,1),2)/norm(ode1_y(:,1),2);
    maxere(iter) = norm(ye1(1,:)'-ode1_y(:,1),inf)/norm(ode1_y(:,1),inf);

    
    rms(iter) = norm(u1(1,:)'-ode1_y(:,1),2)/norm(ode1_y(:,1),2);
    maxer(iter) = norm(u1(1,:)'-ode1_y(:,1),inf)/norm(ode1_y(:,1),inf);
    iter = iter +1
end

figure('name','RMS');
loglog(hspace, rms,'-o');
hold on
loglog(hspace, rmse);
legend('RMS GAUSS-LEGENDRE','RMS EULER')

figure('name','MAX ERROR');
loglog(hspace, maxer,'-o');
hold on
loglog(hspace, maxere);
legend('MAX GAUSS-LEGENDRE','MAX EULER')


function dydt = fun(t,y)
dydt = [y(2);-1.25*y(1) - y(2)];
end