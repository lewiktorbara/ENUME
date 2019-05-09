clear all
close all
[alfa3, alfaValues3, detValues3, condValues3] = minAlfa(10e-5, 10e-5, 3);
[alfa10, alfaValues10, detValues10, condValues10] = minAlfa(10e-5, 10e-5, 10);
[alfa20, alfaValues20, detValues20, condValues20] = minAlfa(10e-5, 10e-5, 20);


%-----


figure('Name', 'Det(A)')
set(gca, 'YScale', 'log')
plot(alfaValues3, detValues3)
hold on

set(gca, 'YScale', 'log')
plot(alfaValues10, detValues10)
hold on

set(gca, 'YScale', 'log')
plot(alfaValues20, detValues20)
hold on
title('Det(A)')
xlabel('alfa')
ylabel('det(A)')
legend('N =3', 'N = 10', 'N = 20')

hold off
%------

figure('Name', 'Cond(A)')
set(gca, 'YScale', 'log')
plot(alfaValues3, condValues3)
hold on

set(gca, 'YScale', 'log')
plot(alfaValues10, condValues10)
hold on

set(gca, 'YScale', 'log')
plot(alfaValues20, condValues20)
hold on
title('Cond(A)')
xlabel('alfa')
ylabel('det(A)')
legend('N =3', 'N = 10', 'N = 20')

hold off
%--------------------------------

cnt=0;
for n=[3,10,20]
    cnt=cnt+1;
    A{1,cnt}=createMatrix(n,1/300);
    B{1,cnt}=invUL(A{1,cnt});
    C{1,cnt}=invLLT(A{1,cnt});
    I=eye(n);
    ARMS(1,cnt) = RMSEr(A{1,cnt}*inv(A{1,cnt}) - I);
    BRMS(1,cnt) = RMSEr(A{1,cnt}*B{1,cnt} - I);
    CRMS(1,cnt) = RMSEr(A{1,cnt}*C{1,cnt} - I);
    AMax(1,cnt) = MaxEr(A{1,cnt}*inv(A{1,cnt}) - I);
    BMax(1,cnt) = MaxEr(A{1,cnt}*B{1,cnt} - I);
    CMax(1,cnt) = MaxEr(A{1,cnt}*C{1,cnt} - I);


    for k=2:22
        A{k,cnt} = A{k-1,cnt};
        A{k,cnt}(:,1) = A{k,cnt}(:,1)*2;
        A{k,cnt}(1,:) = A{k,cnt}(1,:)*2;
        B{k,cnt}=invUL(A{k,cnt});
        C{k,cnt}=invLLT(A{k,cnt});
        ARMS(k,cnt) = RMSEr(A{k,cnt}*inv(A{k,cnt}) - I);
        BRMS(k,cnt) = RMSEr(A{k,cnt}*B{k,cnt} - I);
        CRMS(k,cnt) = RMSEr(A{k,cnt}*C{k,cnt} - I);
        AMax(k,cnt) = MaxEr(A{k,cnt}*inv(A{k,cnt}) - I);
        BMax(k,cnt) = MaxEr(A{k,cnt}*B{k,cnt} - I);
        CMax(k,cnt) = MaxEr(A{k,cnt}*C{k,cnt} - I);
    end
end

for k=0:21
    x(k+1) = 2^k/300;
end

%-------
cnt = 1;
for k=[3,10,20]
    figure('Name', ['RMS(',num2str(k),')'])
    loglog(x,ARMS(:,cnt))
    hold on

    loglog(x,BRMS(:,cnt),'*')
    hold on

    loglog(x,CRMS(:,cnt),'g--')
    hold on
    title(['RMS(',num2str(k),')'])
    xlabel('x')
    ylabel('RMS')

    hold off
    legend('Matlab Function', 'LU Transform', 'LLT Transform')
    
    
    figure('Name', ['Max(',num2str(k),')'])
    loglog(x,AMax(:,cnt))
    hold on

    loglog(x,BMax(:,cnt),'*')
    hold on

    loglog(x,CMax(:,cnt),'g--')
    hold on
    title(['Max(',num2str(k),')'])
    xlabel('x')
    ylabel('MAX')

    hold off
    legend('Matlab Function', 'LU Transform', 'LLT Transform')
    
    cnt=cnt+1;
end


function er = RMSEr(A)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
er = max(sqrt(abs(eig(transpose(A)*A))));
end

function er = MaxEr(A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
er = max(sum(abs(A),2));
end

function [U,L] = LU(A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
U = A;
n = size(A,1);
L = eye(n);
    for i=1:n-1
        for j=i+1:n
                L(j,i)=U(j,i)/U(i,i);
                U(j,:) = U(j,:)-L(j,i)*U(i,:);                
        end
    end

end

function [T,L] = LLT(A)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n = size(A,1);
L = zeros(n);
    for i=1:n
        L(i,i)=sqrt(A(i,i) - sum(L(i,:).^2));
        for j=(i+1):n
            L(j,i) = (A(j,i) - sum(L(i,:).*L(j,:)))/L(i,i);
             
        end
    end
T = transpose(L);

end

function [A] = invUL(B)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[U,L] = LU(B);
n = size(U,1);
I = eye(n);
    for i=1:n
        X(1,i) = I(1,i);
        for j=2:n
            X(j,i) = (I(j,i)-L(j,1:j-1)*X(1:j-1,i));
        end
        A(n,i) = X(n,i)./U(n,n);
        for j=(n-1):-1:1
            A(j,i) = (X(j,i)-U(j,j+1:n)*A(j+1:n,i))./U(j,j);
        end
    end
end

function [A] = invLLT(B)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[T,L] = LLT(B);
n = size(T,1);
I = eye(n);
    for i=1:n
        X(1,i) = I(1,i)./L(1,1);
        for j=2:n
            X(j,i) = (I(j,i)-L(j,1:j-1)*X(1:j-1,i))./L(j,j);
        end
        A(n,i) = X(n,i)/T(n,n);
        for j=(n-1):-1:1
            A(j,i) = (X(j,i)-T(j,j+1:n)*A(j+1:n,i))./T(j,j);
        end
    end
end

function [alfa, alfaValues, detValues, condValues] = minAlfa(precision,difference,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
alfa = 1;
    %while (det(createMatrix(n,sin(alfa - 1))) > precision)
     %   alfa = alfa + difference;
    %end
    %while(det(createMatrix(n,sin(alfa - 1))) > det(createMatrix(n,sin(alfa - 1 + difference))))
     %  alfa = alfa + difference;
    %end
    
    temp = alfa-0.01;
    i=1;
    while temp < alfa+0.01
        alfaValues(i) = temp;
        detValues(i) = det(createMatrix(n,sin(temp-1)));
        condValues(i) = cond(createMatrix(n,sin(temp-1)));
        i = i+1;
        temp = temp + difference;
    end
end

function A = createMatrix(n,x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    for k=n:-1:2
        A(k,k:n)=(k*4)/9;
        A(k:n,k)=(k*4)/9;    
    end
        
    A(1,2:n)=(-2*x)/3;
    A(2:n,1)=(-2*x)/3;
    A(1,1)=x^2;
end


