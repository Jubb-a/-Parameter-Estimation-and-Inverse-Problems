%Math 503B Homework 2 Problem 1 a-d pages 49-51
clear;clf;
%Define Variables
t=[3.49350000000000;4.28530000000000;5.13740000000000;5.81810000000000;6.86320000000000;8.18410000000000];
x=[6;10.1333000000000;14.2667000000000;18.4000000000000;22.5333000000000;26.6667000000000];
sigma=0.1;

%Part a) Find least squares solution for the model t_0 and s_2, plot the
%data, the fitted model and residuals

%fitted model
one=ones(length(x),1);
g=[one, x];
gt=g';
gtg=inv(gt*g);
mApp=gtg*gt*t;
err=sigma*ones(1,length(x));

%Find the least squares solution
v=(1/sigma)*ones(length(x),1);
w=diag(v);
gw=w*g;
mL2=inv(gt*w*w*g)*gt*w*w*t

%plotting data, fitted model and residuals

figure(1);
plot(x,t,'*r'); hold on;
plot(x,mApp(1)+mApp(2)*x,'-ob');
errorbar(x,t,err);
xlabel('x (km)');
ylabel('t (sec)');
legend('Observed Data','Model Data', 'location','best')
title('Part A: Observed Data with Error vs. Model')
%Part b) Calculate and comment on the model parameter correlation matrix

covmL2=inv(transpose(gw)*g)*transpose(gw)*transpose(inv(transpose(gw)*g)*transpose(gw));
corrmatrix=covmL2/sqrt(sigma^2)
%The correlation matrix is positive down the main diagonal and negative in
%the off diagonal. 

%Part c) Plot the error ellipsoid in the (t0,s2) plane. How are the
%correlations manifested in the general appearance of the error ellipsoid
%in (t0,s2) space? Calculate conservative 95% confidence intervals for t0
%and s2 for the appropriate value of chisquared
DELTA2= sum((t-g*mL2).^2/sigma.^2)
plot_ellipse(DELTA2,corrmatrix,mL2)
ci=1.96*sqrt(diag(covmL2))
corrupper=mL2+ci;
corrlower=mL2-ci;
corrint=[corrlower corrupper]
%Since the projetion is narrow with a slight negative slope along its long
%principle axis there is a strong relationship between t_0 and s_2.


%Part d) Evaluate the p-value for this model.
dof=4;
pval=1-chi2cdf(DELTA2,dof)

%PART C FUNCTION CODE
%set the number of points on the ellipse to generate and plot
function plot_ellipse(DELTA2,c,mL2)
n=100;
%construct a vector of n equally-spaced angles from (0,2*pi)
theta=linspace(0,2*pi,n)';
%corresponding unit vector
xhat=[cos(theta),sin(theta)];
cinv=inv(c);
%preallocate output array
r=zeros(n,2);
for i=1:n
%store each (x,y) pair on the confidence ellipse
%in the corresponding row of r
r(i,:)=sqrt(DELTA2/(xhat(i,:)*cinv*xhat(i,:)'))*xhat(i,:);
end

% Plot the ellipse and set the axes.
figure(2);
plot(mL2(1)+r(:,1), mL2(2)+r(:,2));
axis equal
xlabel('t_0 (sec)')
ylabel('s_2(sec/km)')
title('Part C: Ellipsoid Plot')
end

