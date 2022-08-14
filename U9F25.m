%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This program solves Poisson equation
%       u_{xx} + u_{yy} + Ku = f(x,y),  a < x < b,  c < y < d.
%    with Dirichlet boundary condition along  x=a, x=b, y=c, y=d.
%
%    Computed result:% Get 6th order coefficients for interior points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

syms h kw real
assume(h>0)

y=-h:h:h; y=y'; 
x=-h:h:h; x=x';

b=zeros(29,1);
A(1,1) = h;
for i=1:9,
    A(1,i) = 1;
end

for j=-1:1:1,
    for i=-1:1:1,
        k = 1 + i+1 +(j+1)*3;
        xt = x(i+2); 
        yt = y(j+2);
        
        A(2,k) = xt;
        A(3,k) = yt;
        
        A(4,k) = xt*xt/2;
        A(5,k) = yt*yt/2;
        A(6,k) = xt*yt;
        
        A(7,k)  = xt^3/6;
        A(8,k)  = xt^2*yt/2;
        A(9,k)  = xt*yt^2/2;
        A(10,k) = yt^3/6;
        
        A(12,k) = xt^4/24;
        A(13,k) = 4*xt^3*yt/24;
        A(14,k) = 6*xt^2*yt^2/24;
        A(15,k) = 4*xt*yt^3/24;
        A(16,k) = yt^4/24;
        
        A(17,k) = xt^5/120;
        A(18,k) = 5*xt^4*yt/120;
        A(19,k) = 10*xt^3*yt^2/120;
        A(20,k) = 10*xt^2*yt^3/120;
        A(21,k) = 5*xt*yt^4/120;
        A(22,k) = yt^5/120;
        
        A(23,k) = xt^6/720;
        A(24,k) = 6*xt^5*yt/720;
        A(25,k) = 15*xt^4*yt^2/720;
        A(26,k) = 20*xt^3*yt^3/720;
        A(27,k) = 15*xt^2*yt^4/720;
        A(28,k) = 6*xt*yt^5/720;
        A(29,k) = yt^6/720;
    end
end

y=-h:h/2:h; y=y'; 
x=-h:h/2:h; x=x';

for j=-2:1:2,
    for i=-2:1:2,
        k1 = 10 + i+2 +(j+2)*5;
        xt = x(i+3); yt=y(j+3);
        A(1,k1) = - kw;
        
        A(2,k1) = - kw*xt;  
        A(3,k1) = - kw*yt;  
        
        A(4,k1) =  -1 - kw*xt^2/2; 
        A(5,k1) =  -1 - kw*yt^2/2;  
        A(6,k1) =  -kw*xt*yt;  
        
        A(7,k1)  = - xt - kw*xt^3/6;  
        A(8,k1)  = - yt - kw*xt^2*yt/2;  
        A(9,k1)  = - xt - kw*xt*yt^2/2; 
        A(10,k1) = - yt - kw*yt^3/6;  
        
        %beta1 + ...+ beta9 = 1
        A(11,k1) = 1.0;

        A(12,k1) = - (1/2)*xt^2 - kw*xt^4/24;             
        A(13,k1) = - xt*yt - kw*4*xt^3*yt/24;                 
        A(14,k1) = - (1/2)*xt^2 - (1/2)*yt^2 - kw*6*xt^2*yt^2/24;    
        A(15,k1) = - xt*yt - kw*4*xt*yt^3/24;                   
        A(16,k1) = - (1/2)*yt^2 - kw*yt^4/24;             
        
        A(17,k1) = - (1/6)*xt^3 - kw*xt^5/120;                 
        A(18,k1) = - (1/2)*xt^2*yt - kw*5*xt^4*yt/120 ;               
        A(19,k1) = - (1/6)*xt^3 - (1/2)*xt*yt^2 - kw*10*xt^3*yt^2/120;     
        A(20,k1) = - (1/6)*yt^3 - (1/2)*xt^2*yt - kw*10*xt^2*yt^3/120 ;    
        A(21,k1) = - (1/2)*xt*yt^2 - kw*5*xt*yt^4/120 ;              
        A(22,k1) = - (1/6)*yt^3 - kw*yt^5/120;                 

        A(23,k1) = - (1/24)*xt^4;                    
        A(24,k1) = - (1/6)*xt^3*yt ;                  
        A(25,k1) = - (1/24)*xt^4 - (1/4)*xt^2*yt^2;      
        A(26,k1) = - (1/6)*xt^3*yt - (1/6)*xt*yt^3;     
        A(27,k1) = - (1/4)*xt^2*yt^2 - (1/24)*yt^4;      
        A(28,k1) = - (1/6)*xt*yt^3;                  
        A(29,k1) = - (1/24)*yt^4;                    
    end
end

b(11) = 1;
Ao = A; bo = b;

ind1 = 9 + [2,4,6,8,10,12,14,16,18,20,22,24];
A(:,ind1) = [];

Rank = [rank(A), rank([A,b])]
xs = A\b;

xs1 = zeros(34,1)*h;
ind3 = setdiff(1:34,ind1);
xs1(ind3(:)) = xs;

UC = reshape(xs1(1:9)', 3, 3)
FC = reshape(xs1(10:end)', 5, 5)
