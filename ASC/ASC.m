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

  
opts = odeset('Reltol',exp(-13),'AbsTol',exp(-14),'Stats' ,'on');
t = linspace(0,10,1001);
[ode_t, ode_y] = ode113(@fun,t,y0',opts);


  
iterator = 2;
for i=0.01:h:10
    F(:,iterator) = Bu\[A*u(:,iterator-1);A*u(:,iterator-1) ;A*u(:,iterator-1)];
    u(:,iterator) = u(:,iterator-1) + h*((5/18)*F(1:2,iterator) + (4/9)*F(3:4,iterator) +(5/18)*F(5:6,iterator));
        
    ye(:,iterator) = ye(:,iterator-1) + A*h*ye(:,iterator-1);
    
    iterator = iterator +1;
    
end

figure('name','1')
plot(t,u(1,:));
hold on
plot(t, ode_y(:,1));
hold off

figure('name','2')
plot(t,ye(1,:));
hold on
plot(t, ode_y(:,1));
hold off



%------------------------------------------------------------   
hspace = logspace(-4,-1,100);
iter = 1;    
for h=hspace
    clear u1 ode1_y t1 ode1_t F ye1;
    
    u1(:,1) = [4 -2];
    ye1(:,1) = [4 -2];
    Bu = [I-(5/36)*h*A, ((sqrt(15)/15)-(2/9))*h*A, ((sqrt(15)/30)-5/36)*h*A;
      -(5/36 + sqrt(15)/24)*h*A, I-(2/9)*h*A, ((sqrt(15)/24)-5/36)*h*A; 
      -(5/36 + sqrt(15)/30)*h*A, -(2/9 + sqrt(15)/15)*h*A, I-(5/36)*h*A];


    t1 = linspace(0,10,size(0.01:h:10,2)+1);
    [ode1_t, ode1_y] = ode113(@fun,t1,y0',opts);


  
    iterator = 2;
    for i=0.01:h:10
        F(:,iterator) = Bu\[A*u1(:,iterator-1);A*u1(:,iterator-1) ;A*u1(:,iterator-1)];
        u1(:,iterator) = u1(:,iterator-1) + h*((5/18)*F(1:2,iterator) + (4/9)*F(3:4,iterator) +(5/18)*F(5:6,iterator));
        
        ye1(:,iterator) = ye1(:,iterator-1) + A*h*ye1(:,iterator-1);
        
        iterator = iterator +1;
        
    end
      
    rmse(iter) = norm(ye1(1,:)'-ode1_y(:,1),2)/norm(ode1_y(:,1),2);
    maxere(iter) = norm(ye1(1,:)'-ode1_y(:,1),inf)/norm(ode1_y(:,1),inf);

    
    rms(iter) = norm(u1(1,:)'-ode1_y(:,1),2)/norm(ode1_y(:,1),2);
    maxer(iter) = norm(u1(1,:)'-ode1_y(:,1),inf)/norm(ode1_y(:,1),inf);
    iter = iter +1
end

figure('name','RMS');
loglog(hspace, rms);

figure('name','MAX ERROR');
loglog(hspace, maxer);

figure('name','RMS EULER');
loglog(hspace, rmse);

figure('name','MAX ERROR EULER');
loglog(hspace, maxere);


function dydt = fun(t,y)
dydt = [y(2);-1.25*y(1) - y(2)];
end