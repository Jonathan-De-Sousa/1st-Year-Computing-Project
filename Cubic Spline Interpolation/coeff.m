function polycoeff=coeff(x,f)
%coeff - Function used to evaluate the polynomial coefficients of an interpolant function
%for each interval of a data set consisting of arrays x and f
%INPUTS 
%x - Vector containing independent variable values for data set
%f - Vector of corresponding dependent variable values for data set
%OUTPUTS
%polycoeff - P x 4 matrix containing the coefficients of the interpolant 
%            function where P is the number of intervals.
    
%Calculate number of x values
N=length(x);

%Initialise matrices A, h, b and polycoeff
%NOTE: the problem will be set up as a system Am=b 
%where m - matrix containing curvature values of the interpolant function 
%          at each each x value
%      h - matrix containing width of intervals
%      A - matrix containing combinations of h values
%      b - matrix containing combinations of y values
h=zeros(1,(N-1));     
A=zeros((N-2),N);     
b=zeros((N-2),1);       
polycoeff=zeros((N-1),4);

%Calculate width of the each interval
for i=1:N-1
    h(i)= x(i+1)-x(i);
end

%Calculate elements of matrix A
for i=2:N-1
    A((i-1),(i-1))=h(i-1);
    A((i-1),i)=2*(h(i-1)+h(i));
    A((i-1),(i+1))=h(i);
end

%Calculate elements of matrix b
for i=2:N-1
    b(i-1)=6*((f(i+1)-f(i))/h(i)-((f(i)-f(i-1))/h(i-1)));
end

%Since the boundary conditions for m, the interpolant function curvature,
%is m(1)=0 and m(N)=0, the system of equations can be truncated
%Remove first and last columns of matrix A
A=A(:,2:(N-1));  

%Calculate matrix m 
m=A\b;

%Add the external boundary condition values of m
m=[0;m;0];

%Calculate polycoeff
for i=1:N-1
    polycoeff(i,4)=f(i);    %constant term coefficient
    polycoeff(i,3)=((f(i+1)-f(i))/h(i))-(h(i)*(2*m(i)+m(i+1))/6);   %linear term coefficient
    polycoeff(i,2)=m(i)/2;  %quadratic term coefficient
    polycoeff(i,1)=(m(i+1)-m(i))/(6*h(i));  %cubic term coefficient
end

end