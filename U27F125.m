%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This program solves Poisson equation
%       u_{xx} + u_{yy} + u_{zz} + Ku = f(x,y,z),
%    with Dirichlet boundary condition along x=a, x=b, y=c, y=d.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

syms h sigma kw
assume(h>0)
assume(sigma>0)

y=-h:h:h; y=y'; x=-h:h:h; x=x'; z=-h:h:h; z=z';
b=zeros(85,1);

A(1,1) = h;
for i=1:27
    A(1,i)=1;
end

for k=-1:1:1,  
    for j=-1:1:1,
        for i=-1:1:1,  
            k1 = 9*(k + 2 - 1) + 3*(j + 2 - 1) + (i + 2);
            xt = x(i+2); yt = y(j+2); zt = z(k+2);

            A(2,k1) = xt; A(3,k1) = yt; A(4,k1) = zt;
            
            A(5,k1) = xt^2/2;  A(6,k1) = yt^2/2;     A(7,k1) = zt^2/2;
            A(8,k1) = xt*yt;    A(9,k1) = xt*zt;     A(10,k1) = yt*zt;
            
            A(11,k1) = xt^3/6; 	A(12,k1) = yt^3/6; 	A(13,k1) = zt^3/6;
            A(14,k1) = xt*xt*yt/2; A(15,k1) = xt*yt*yt/2; A(16,k1) = xt*xt*zt/2; 
            A(17,k1) = xt*zt*zt/2; A(18,k1) = yt*yt*zt/2; A(19,k1) = yt*zt*zt/2;                      
            A(20,k1) = xt*yt*zt;      
            
            A(21,k1) = xt^4/24;          A(22,k1) = yt^4/24;          A(23,k1) = zt^4/24;
            A(24,k1) = 4*xt^3*yt/24;     A(25,k1) = 4*xt*yt^3/24;     A(26,k1) = 4*xt^3*zt/24;        
            A(27,k1) = 4*xt*zt^3/24;     A(28,k1) = 4*yt^3*zt/24;     A(29,k1) = 4*yt*zt^3/24;
            A(30,k1) = 6*xt^2*yt^2/24;   A(31,k1) = 6*xt^2*zt^2/24;   A(32,k1) = 6*yt^2*zt^2/24;      
            A(33,k1) = 12*xt^2*yt*zt/24; A(34,k1) = 12*xt*yt^2*zt/24; A(35,k1) = 12*xt*yt*zt^2/24;    
            
            A(37,k1) = xt^5/120;            A(38,k1) = yt^5/120;            A(39,k1) = zt^5/120;
            A(40,k1) = 5*xt^4*yt/120;       A(41,k1) = 5*xt*yt^4/120;       A(42,k1) = 5*xt^4*zt/120;
            A(43,k1) = 5*xt*zt^4/120;       A(44,k1) = 5*yt^4*zt/120;       A(45,k1) = 5*yt*zt^4/120; 
            A(46,k1) = 10*xt^3*yt^2/120;    A(47,k1) = 10*xt^2*yt^3/120;    A(48,k1) = 10*xt^3*zt^2/120;
            A(49,k1) = 10*xt^2*zt^3/120;    A(50,k1) = 10*yt^3*zt^2/120;    A(51,k1) = 10*yt^2*zt^3/120;
            A(52,k1) = 20*xt^3*yt*zt/120;   A(53,k1) = 20*xt*yt^3*zt/120;   A(54,k1) = 20*xt*yt*zt^3/120;
            A(55,k1) = 30*xt^2*yt^2*zt/120; A(56,k1) = 30*xt^2*yt*zt^2/120; A(57,k1) = 30*xt*yt^2*zt^2/120;
            
            A(58,k1) = xt^6/720;           A(59,k1) = yt^6/720;         A(60,k1) = zt^6/720;         
            A(61,k1) = xt^5*yt/120;        A(62,k1) = yt^5*xt/120;      A(63,k1) = xt^5*zt/120;
            A(64,k1) = xt*zt^5/120;        A(65,k1) = yt^5*zt/120;      A(66,k1) = yt*zt^5/120;       
            A(67,k1) = xt^4*yt^2/48;       A(68,k1) = xt^2*yt^4/48;     A(69,k1) = xt^4*zt^2/48;
            A(70,k1) = xt^2*zt^4/48;       A(71,k1) = yt^4*zt^2/48;     A(72,k1) = yt^2*zt^4/48;  
            A(73,k1) = xt^3*yt^3/36;       A(74,k1) = xt^3*zt^3/36;     A(75,k1) = yt^3*zt^3/36;     
            A(76,k1) = xt^4*yt*zt/24;      A(77,k1) = xt*yt^4*zt/24;    A(78,k1) = xt*yt*zt^4/24; 
            A(79,k1) = xt^3*yt^2*zt/12;    A(80,k1) = xt^2*yt^3*zt/12;  A(81,k1) = xt^3*yt*zt^2/12;
            A(82,k1) = xt^2*yt*zt^3/12;    A(83,k1) = xt*yt^3*zt^2/12;  A(84,k1) = xt*yt^2*zt^3/12;
            A(85,k1) = xt^2*yt^2*zt^2/8;
        end
    end
end

x=-h:h/2:h; x=x';
y=-h:h/2:h; y=y'; 
z=-h:h/2:h; z=z'; 
for  k=-2:1:2,
    for j=-2:1:2,
        for i=-2:1:2,
            k2 = 27 + 25*(k+2) + 5*(j+2) + (i+3);
            xt = x(i+3); yt=y(j+3); zt=z(k+3);       

            A(1,k2) = 0 - kw;   
            
            A(2,k2) = 0 - kw*xt;   
            A(3,k2) = 0 - kw*yt;  
            A(4,k2) = 0 - kw*zt;
            
            A(5,k2) = -1 - kw*xt^2/2;   
            A(6,k2) = -1 - kw*yt^2/2;  
            A(7,k2) = -1 - kw*zt^2/2;         
            A(8,k2) =  0 - kw*xt*yt;
            A(9,k2) =  0 - kw*xt*zt;
            A(10,k2) = 0 - kw*yt*zt;

            A(11,k2) = -xt - kw*xt^3/6;
            A(12,k2) = -yt - kw*yt^3/6;
            A(13,k2) = -zt - kw*zt^3/6;
            A(14,k2) = -yt - kw*xt^2*yt/2;
            A(15,k2) = -xt - kw*xt*yt^2/2;
            A(16,k2) = -zt - kw*xt^2*zt/2;
            A(17,k2) = -xt - kw*xt*zt^2/2;
            A(18,k2) = -zt - kw*yt^2*zt/2;    
            A(19,k2) = -yt - kw*yt*zt^2/2;    
            A(20,k2) = -kw*xt*yt*zt;      

            A(21,k2) = A(21,k2) - xt^2/2 - kw*xt^4/24;
            A(22,k2) = A(22,k2) - yt^2/2 - kw*yt^4/24;
            A(23,k2) = A(23,k2) - zt^2/2 - kw*zt^4/24;
            A(24,k2) = A(24,k2) - xt*yt - kw*xt^3*yt/6;
            A(25,k2) = A(25,k2) - xt*yt - kw*xt*yt^3/6;
            A(26,k2) = A(26,k2) - xt*zt - kw*xt^3*zt/6;
            A(27,k2) = A(27,k2) - xt*zt - kw*xt*zt^3/6;
            A(28,k2) = A(28,k2) - yt*zt - kw*yt^3*zt/6;
            A(29,k2) = A(29,k2) - yt*zt - kw*yt*zt^3/6;
            A(30,k2) = A(30,k2) - xt^2/2 - kw*xt^2*yt^2/4;  A(30,k2) = A(30,k2) - yt^2/2;                                       
            A(31,k2) = A(31,k2) - xt^2/2 - kw*xt^2*zt^2/4;  A(31,k2) = A(31,k2) - zt^2/2;            
            A(32,k2) = A(32,k2) - yt^2/2 - kw*yt^2*zt^2/4;  A(32,k2) = A(32,k2) - zt^2/2;
            A(33,k2) = A(33,k2) - yt*zt - kw*xt^2*yt*zt/2;
            A(34,k2) = A(34,k2) - xt*zt - kw*xt*yt^2*zt/2;       
            A(35,k2) = A(35,k2) - xt*yt - kw*xt*yt*zt^2/2;        
            A(36,k2) = 1.0;                     
                                                
            A(37,k2) = A(37,k2)-xt^3/6 - kw*xt^5/120;
            A(38,k2) = A(38,k2)-yt^3/6 - kw*yt^5/120;
            A(39,k2) = A(39,k2)-zt^3/6 - kw*zt^5/120;
            
            A(40,k2) = A(40,k2)-xt^2*yt/2 - kw*xt^4*yt/24;          
            A(41,k2) = A(41,k2)-yt^2*xt/2 - kw*xt*yt^4/24;
            A(42,k2) = A(42,k2)-xt^2*zt/2 - kw*xt^4*zt/24;
            A(43,k2) = A(43,k2)-xt*zt^2/2 - kw*xt*zt^4/24;
            A(44,k2) = A(44,k2)-yt^2*zt/2 - kw*yt^4*zt/24;
            A(45,k2) = A(45,k2)-yt*zt^2/2 - kw*yt*zt^4/24;
            
            A(46,k2) = A(46,k2)-xt^3/6;     A(46,k2) = A(46,k2)-yt^2*xt/2 - kw*xt^3*yt^2/12;
            A(47,k2) = A(47,k2)-yt^3/6;     A(47,k2) = A(47,k2)-xt^2*yt/2 - kw*xt^2*yt^3/12;
            A(48,k2) = A(48,k2)-xt^3/6;     A(48,k2) = A(48,k2)-xt*zt^2/2 - kw*xt^3*zt^2/12;
            A(49,k2) = A(49,k2)-zt^3/6;     A(49,k2) = A(49,k2)-xt^2*zt/2 - kw*xt^2*yt^3/12;
            A(50,k2) = A(50,k2)-yt^3/6;     A(50,k2) = A(50,k2)-yt*zt^2/2 - kw*yt^3*zt^2/12;
            A(51,k2) = A(51,k2)-zt^3/6;     A(51,k2) = A(51,k2)-yt^2*zt/2 - kw*yt^2*zt^3/12;
            A(52,k2) = A(52,k2)-xt*yt*zt - kw*xt^3*yt*zt/6;
            A(53,k2) = A(53,k2)-xt*yt*zt - kw*xt*yt^3*zt/6;
            A(54,k2) = A(54,k2)-xt*yt*zt - kw*xt*yt*zt^3/6;
            A(55,k2) = A(55,k2)-xt^2*zt/2;  A(55,k2) = A(55,k2)-yt^2*zt/2 - kw*xt^2*yt^2*zt/4;
            A(56,k2) = A(56,k2)-xt^2*yt/2;  A(56,k2) = A(56,k2)-yt*zt^2/2 - kw*xt^2*yt*zt^2/4;
            A(57,k2) = A(57,k2)-yt^2*xt/2;  A(57,k2) = A(57,k2)-xt*zt^2/2 - kw*xt*yt^2*zt^2/4;

            A(58,k2) = A(58,k2)-xt^4/24;
            A(59,k2) = A(59,k2)-yt^4/24;
            A(60,k2) = A(60,k2)-zt^4/24;
            
            A(61,k2) = A(61,k2)-xt^3*yt/6;
            A(62,k2) = A(62,k2)-xt*yt^3/6;
            A(63,k2) = A(63,k2)-xt^3*zt/6;
            A(64,k2) = A(64,k2)-xt*zt^3/6;
            A(65,k2) = A(65,k2)-yt^3*zt/6;
            A(66,k2) = A(66,k2)-yt*zt^3/6;
            A(67,k2) = A(67,k2)-xt^4/24;     A(67,k2) = A(67,k2)-xt^2*yt^2/4;
            A(68,k2) = A(68,k2)-yt^4/24;     A(68,k2) = A(68,k2)-xt^2*yt^2/4;
            A(69,k2) = A(69,k2)-xt^4/24;     A(69,k2) = A(69,k2)-xt^2*zt^2/4;
            
            A(70,k2) = A(70,k2)-zt^4/24;     A(70,k2) = A(70,k2)-xt^2*zt^2/4;
            A(71,k2) = A(71,k2)-yt^4/24;     A(71,k2) = A(71,k2)-yt^2*zt^2/4;
            A(72,k2) = A(72,k2)-zt^4/24;     A(72,k2) = A(72,k2)-yt^2*zt^2/4;
            A(73,k2) = A(73,k2)-xt^3*yt/6;   A(73,k2) = A(73,k2)-xt*yt^3/6;
            
            A(74,k2) = A(74,k2)-xt*zt^3/6;   A(74,k2) = A(74,k2)-xt^3*zt/6;
            A(75,k2) = A(75,k2)-yt^3*zt/6;   A(75,k2) = A(75,k2)-yt*zt^3/6;
            
            A(76,k2) = A(76,k2)-xt^2*yt*zt/2;
            A(77,k2) = A(77,k2)-xt*yt^2*zt/2;
            A(78,k2) = A(78,k2)-xt*yt*zt^2/2;
            A(79,k2) = A(79,k2)-xt^3*zt/6;    A(79,k2) = A(79,k2)-xt*yt^2*zt/2;
            A(80,k2) = A(80,k2)-yt^3*zt/6;    A(80,k2) = A(80,k2)-xt^2*yt*zt/2;
            A(81,k2) = A(81,k2)-xt^3*yt/6;    A(81,k2) = A(81,k2)-xt*yt*zt^2/2;
            A(82,k2) = A(82,k2)-yt*zt^3/6;    A(82,k2) = A(82,k2)-xt^2*yt*zt/2;
            A(83,k2) = A(83,k2)-xt*yt^3/6;    A(83,k2) = A(83,k2)-xt*yt*zt^2/2;
            
            A(84,k2) = A(84,k2)-xt*zt^3/6;    A(84,k2) = A(84,k2)-xt*yt^2*zt/2;
            A(85,k2) = A(85,k2)-yt^2*zt^2/4;  A(85,k2) = A(85,k2)-xt^2*zt^2/4;   A(85,k2) = A(85,k2)-xt^2*yt^2/4;
        end
    end
end
Ao = A; b(36) = 1;

ind1 = 27 +       [1,2,4,5,6,7,8,9,10,12,14,16,17,18,19,20,21,22,24,25];
ind2 = 27 + 25 +  [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25];
ind3 = 27 + 50 +  [2,4,6,7,9,10,16,17,19,20,22,24];
ind4 = 27 + 75 +  [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25];
ind5 = 27 + 100 + [1,2,4,5,6,7,8,9,10,12,14,16,17,18,19,20,21,22,24,25];
ind = [ind1,ind2,ind3,ind4,ind5];
A(:,ind) = [];

Rank = [rank(A), rank([A,b])];

xs = Calculate(A,b);

xs1 = zeros(152,1)*h;
indn = setdiff(1:152,ind);
xs1(indn,:) = xs;

UC = xs1(1:27)'
FC = xs1(28:end)'