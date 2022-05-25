%Computation task 3. relationship between variance and sales-weighted
%markup

%variance ={1, 4, 9, 16, 25}.
for A = [1,2,3,4,5]  
i(A)=A;
I=3
J=100;%should be 100 
P=1;
L=1;
t=1.5; %theta
n=10;%ita
a=lognrnd(0,A, [J I]);

%create marketshare matrix
x=sym('x', [J I]);assumeAlso(x < 1 & x>0);assumeAlso(sum(x,2)==1);

%equation based on computation task 1.b
M = ((1./(a .*(1- 1/n +(1/n - 1/t).*x))).^(1-n));z=sum(M,2);c=M./z;

fun1=matlabFunction(x-c, 'vars', {x} );
x0=ones(J,I);
options = optimoptions('fsolve','Display','none');
s=fsolve(fun1,x0,options)

mu = 1./(1- 1/n + (1/n - 1/t).*s);%markup

%wage
W1= (sum(sum( (1/I)*(a/mu).^(n-1),2).^((t-1)/(n-1))*(1/J),"all"))^(1/(t-1))

mc=W1./a; %marginal cost
p=mc.*mu; %price

%create firm output matrix
yij=sym('yij', [J I]);
assume(yij>0);
assumeAlso(sum(yij./a,2)==1)

%based on the marketshare defination 
TP=sum(p.*yij,2)
output=s.*TP./p

fun2=matlabFunction(output-yij, 'vars', {yij} );
y0=ones(J,I);
options = optimoptions('fsolve','Display','none');
y=fsolve(fun1,y0,options)

Yj= sum(I^(-1/n).*y.^((n-1)/n),2).^(n/(n-1)) %Yj market level output
Y=sum(J^(-1/t).*Yj.^((t-1)/t)).^(t/(t-1)); %Y aggregate output

a_MU(A)=sum((p.*y)./(P.*Y).*mu,"all")/J;%sales weighted average markup

 end

 
 
 plot(i,a_MU)
 xlabel('I')
ylabel('sales weighted average markups')


