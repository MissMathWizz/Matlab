%Computation task 1 and 2 
for I = [2,3,5,7,10]  
i(I)=I;
J=100;%should be 100 
P=1;
L=1;
t=1.5; %theta
n=10;%ita
a=lognrnd(0,2, [J I]);

%create marketshare s matrix
x=sym('x', [J I]);assumeAlso(x < 1 & x>0);assumeAlso(sum(x,2)==1);

%replicates the equation in computation task 1.b
M = ((1./(a .*(1- 1/n +(1/n - 1/t).*x))).^(1-n));z=sum(M,2);c=M./z;

fun1=matlabFunction(x-c, 'vars', {x} );
x0=ones(J,I);
options = optimoptions('fsolve','Display','none');
s=fsolve(fun1,x0,options)


mu = 1./(1- 1/n + (1/n - 1/t).*s);%markup

%Aggregate Wage based on 2021 paper appendix
W1= (sum(sum( (1/I)*(a/mu).^(n-1),2).^((t-1)/(n-1))*(1/J),"all"))^(1/(t-1))
vpa(W1)

mc=W1./a;%marginal cost
p=mc.*mu;%price

%create price matrix
yij=sym('yij', [J I]);
assume(yij>0);
assumeAlso(sum(yij./a,2)==1)

%Based on marketshare defination
TP=sum(p.*yij,2);%total profit
output=s.*TP./p;% marketshare = price * Y/ total profit

fun2=matlabFunction(output-yij, 'vars', {yij} );
y0=ones(J,I);
options = optimoptions('fsolve','Display','none');
y=fsolve(fun1,y0,options)

Yj= sum(I^(-1/n).*y.^((n-1)/n),2).^(n/(n-1)) %Yj market level output
Y(I)=sum(J^(-1/t).*Yj.^((t-1)/t)).^(t/(t-1)); %Y aggregate output

a_MU(I)=sum((p.*y)/(P.*Y).*mu,"all")/J;%sales weighted average markup
Profit = sum(p.*y-(W1./a).*y,"all");
Profit_j(I)= Profit/J; %economy-wide average profit
Wage(I)=W1 %wage;
max_a_Profit_I(I)= Profit/J/I;

 end



