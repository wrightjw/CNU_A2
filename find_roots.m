clear all,close all %Clearing all prior variables, Closing all prior graphs
optimoptions('fsolve','TolFun',1e-10); %Changing the tolerance on fsolve - taken from https://www.mathworks.com/help/optim/ug/fsolve.html
format short %Formatiing output where not specified

%3.1
x_plot = linspace(-1,1,1000); %Generating 1000 evenly spaced values between -1 and 1 to plot
F_plot = (cos(pi.*x_plot)).^3+0.01.*(x_plot+2).^(1/10); %Generating y-axis values to plot
figure(1); %Creating first figure
plot(x_plot,F_plot,'-r'); %Plotting values and specifying line type
legend('F(x)'); %Creating legend
title('Graph of $F(x)=\cos^{3}(\pi x)+0.01(x+2)^{\frac{1}{10}}$ against $x$','interpreter','latex'); %Titling graph
xlim([min(x_plot) max(x_plot)]); %Setting x axis limits 
ylim([min(F_plot) max(F_plot)]); %Setting y axis limits
xlabel('$x$','interpreter','latex'); %Labelling x axis
ylabel('$F(x)$','interpreter','latex'); %Labelling y axis
F = @(x) cos(pi*x)^3+0.01*(x+2)^(1/10); %Defining function F from (6) for fsolve to use
x_star_negative = fsolve(F,-0.6,optimset('display','off')); %Calculating the negative root by starting at -0.6, guessed from graph, and turning output of fsolve off - optimset taken from https://uk.mathworks.com/matlabcentral/answers/2583-how-to-remove-text-from-fsolve-result
x_star_positive = fsolve(F,0.4,optimset('display','off')); %Calculating the positive root by starting at 0.4, guessed from graph, and turning output of fsolve off
fprintf('3.1) The negative root is x_star_negative = %.9f and the positive root is x_star_positive = %.9f\n' ,x_star_negative,x_star_positive) %Displaying roots in one line with explantory text to 9 decimal places

%3.2
x_val = -2/3; %Assigning initial value of x
n = 0; %Assigning initial step
F_val = (cos(pi*x_val))^3+0.01*(x_val+2)^(1/10); %Testing whether function is not already s.t. abs(F)<10^(-10) at initial point

while abs(F_val) >= 10^(-10) %Will loop until |F|<10^(-10)
    F_val = (cos(pi*x_val))^3+0.01*(x_val+2)^(1/10); %Calculating function value for the prior x value in this iteration
    x_val = x_val + x_val/abs(x_val) * F_val; %Assigning new value of x using equation (7)
    n = n + 1; %Defining number of iterations taken to find root
end

fprintf('3.2) Number of iterations taken is %u and the negative root value is x = %.9f to an error of 10^(-10). \n',n,x_val) %Displaying number of iterations and the negative  root value to 9 decimal places

%3.3
syms x; %defining variable x as a symbol for differentiation in-built function - looked up at https://uk.mathworks.com/help/symbolic/diff.html
F = (cos(pi*x))^3+0.01*(x+2)^(1/10); %defining function to be differentiated
deriv_F = diff(F,1); %differentiating function using in-built - looked up at https://uk.mathworks.com/help/symbolic/diff.html
x_val = -2/3; %Defining initial x value
n = 0; %Setting initial step
F = (cos(pi*x_val))^3+0.01*(x_val+2)^(1/10); %Evaluating F for x=-2/3 to check not already s.t. abs(F)<10^(-10)

while abs(F) >= 10^(-10) %Looping until |F|<10^(-10)
    F = (cos(pi*x_val))^3+0.01*(x_val+2)^(1/10); %Evaluating F(x) at each x_val calculated
    deriv_val = double(vpa(subs(deriv_F,x,x_val))); %Evaluating exact value of deriv_F at each x_val and turning to double to deal with data-type issue - looked up at https://uk.mathworks.com/help/symbolic/differentiation.html
    x_val = x_val - F/deriv_val; %Assigning new value of x_val using Newton's method
    n = n+1; %Defining number of iterations taken to find root
end

fprintf('3.3) Number of iterations taken is %u and the negative root value is %.9f to an error of 10^(-10). \n',n,x_val) %Displaying output for number of iterations and root to 9 decimal places

%3.4
x_fp = [2/3]; %Vector containing initial value 2/3 for fixed point iteration
x_newton = [2/3]; %Vector containing initial value 2/3 for newton method

for i = 2:50 %Looping values of x_fp and x_newton into respective sets for 100 values (used in 3.5)
    F_fp = (cos(pi*x_fp(i-1)))^3+0.01*(x_fp(i-1)+2)^(1/10); %Calculating function values for fp method
    F_newton = (cos(pi*x_newton(i-1)))^3+0.01*(x_newton(i-1)+2)^(1/10); %Calculating function values for newton method
    x_fp(i) = x_fp(i-1) + x_fp(i-1)/abs(x_fp(i-1))*F_fp; %Calculating next x root using fp method
    deriv_val = double(vpa(subs(deriv_F,x,x_newton(i-1)))); %Reusing derivative calculated in 3.3
    x_newton(i) = x_newton(i-1) - F_newton/deriv_val; %Calculating next x root using newtons method
end

difference_fp = abs(x_fp-x_star_positive); %Magnitude of difference between x_fp and x_star_positive
difference_newton = abs(x_newton-x_star_positive); %Magnitude of difference between x_newton and x_star_positive

figure(2); %Creating second figure
n_plot = linspace(0,6,7); %Creating x axis plot of n=0,1,2,3,4,5,6.
plot(n_plot,difference_fp(1:7),'-g'), hold on, plot(n_plot, difference_newton(1:7),'-b'); %Plotting n vs the difference in x_fp and x_star_positive and plotting n vs the difference in x_newton and x_star_positive corresponding to the n=0,...,6 values
legend('Fixed Point Method Error','Newton Method Error'); %Creating legend
title('Magnitude of Difference Between Newton Method and Root and Fixed Point Method and Root Respectively'); %Titling graph
xlim([min(n_plot) max(n_plot)]); %Setting x axis limits 
ylim([0 max(difference_fp)]); %Setting y axis limits
xlabel('Number of Steps'); %Labelling x-axis
ylabel('Magnitude of Difference'); %Labelling y-axis

%3.5

log_diff_fp = log(difference_fp); %log of difference in fixed point method
log_diff_newton = log(difference_newton); %log of difference in newton method
steps = linspace(0,length(x_newton)-1,length(x_newton)); %Creating x-axis points to plot graph

figure(3); %Creating third figure
plot(steps,log_diff_fp,'--g'), hold on, plot(steps,log_diff_newton,'--r'); %Plotting the log of the differences to show the rate of convergence
legend('Logarithm of Fixed Point Method Magnitude Error', 'Logarithm of Newton Method Magnitude Error'); %Adding legend
title('Logarithm of Magnitude Errors of Methods vs Steps'); %Titling figure
xlim([min(steps) max(steps)]); %Setting x axis limits 
ylim([min(log_diff_newton)-0.5 0]); %Setting y axis limits
xlabel('Number of Steps'); %Labelling x axis
ylabel('Logarithm of Error'); %Labelling y axis

q1 = zeros(1,length(x_fp)-1); %Creating array to loop quotients for fixed point method into
x_plot_quotient = steps(1:length(x_fp)-1); %Creating array for x-axis plot for figure 2
for i=2:length(x_fp)-1 %Looping over differences in newton method x values and fixed point method x values
    q1(i)=log(abs(x_fp(i+1)-x_fp(i)))/log(abs(x_fp(i)-x_fp(i-1))); %Appending quotients for fixed point method - formula from https://en.wikipedia.org/wiki/Rate_of_convergence
end

figure(4); %Creating fourth figure
plot(x_plot_quotient,q1,'g'); %Plotting the convergence quotient
legend('Convergence Quotient - Fixed Point'); %Adding legend
title('Convergence Quotient Plots'); %Titling figure
xlim([min(x_plot_quotient) max(x_plot_quotient)]); %Setting x axis limits 
ylim([0 3]); %Setting y axis limits
xlabel('Step Difference'); %Labelling x axis
ylabel('Convergence Quotient'); %Labelling y axis

format long %Changing number of decimals displayed
disp('3.5) The 3rd to 6th values of the Newton Method approximation magnitude errors:') %Printing explanatory text
Errors_N = difference_newton(3:6) %Displaying the magnitude errors for Newton Method
disp('3.5) The 3rd to 6th values of the Fixed Point Method approximation magnitude errors:') %Printing explanatory text
Errors_FP = difference_fp(3:6) %Displaying the magnitude errors for Newton Method

%Figure(2) illustrates that the Newton Method converges faster than the
%fixed point method. The reason for this is demonstrated by figure(3) whereby the
%gradient for the newton method is greater than the gradient
%for the fixed point method beyond n=2, illustrating the greater rate of convergence.
%Note that the 'Errors_N' output in the command window are approximately
%10^(-2),10^(-3),10^(-5),10^(-6) and so approximately a factor of 10 each time, 
%so the newton method is converging logarithmically.
%Figure(4) illustrates using the rate of convergence approximation formula
%that the Fixed Point method converges logarithmically as the convergence quotient
%is approximately 1 for large n. This plot did not work for Newton Method has
%as the difference between the estimate and true root is too close to zero
%after n>7. We see for 'Errors_FP' that they are decreasing by approximately a factor of 2 each
%time, explaining the difference in gradients in figure(3).
%Further, from the displayed values of 3.2 and 3.3 show us that to converge 
%to the same root approximation within 10^(-10) of accuracy, the fixed 
%point method takes 35 iterations whilst the newton method takes 7 iterations, 
%meaning that the newton method converges on that root 5 times faster in this case.