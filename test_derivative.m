clear all %Clearing all prior variables
close all %Closing all prior graphs
format short %Formatting decimal places

%1.1
x_input = pi/2; %Defining value of x being evaluated
syms x; %Defining variable x as a symbol for differentiation in-built function - looked up at https://uk.mathworks.com/help/symbolic/diff.html
F = exp(x)*sin(x); %Defining function to be differentiated
deriv_F = diff(F,1); %Differentiating function using in-built - looked up at https://uk.mathworks.com/help/symbolic/diff.html
deriv_F = double(vpa(subs(deriv_F,x,x_input))); %Evaluating exact value of deriv_F at input value x_input and turning to double to deal with data-type issue for line 17 - looked up at https://uk.mathworks.com/help/symbolic/differentiation.html

dx = [0.0032,0.0016,0.0008,0.0004,0.0002,0.0001]; %Creating array of dx to evaluate
Error_Magnitude=zeros(1,length(dx)); %Preallocating array to loop error magnitudes into for speed
for i=1:length(dx) %Running over number of values dx in question
    D = (1/dx(i))*((1/28)*(exp(x_input+5*dx(i))*sin(x_input+5*dx(i)))+(1/4)*(exp(x_input+dx(i))*sin(x_input+dx(i)))-(2/7)*((exp(x_input-2*dx(i))*sin(x_input-2*dx(i))))); %Evaluating approximation for each dx, with F(x) substituted using given formula
    Error_Magnitude(i) = abs(deriv_F-D); %Calculating error magnitude
    fprintf(' %u) Approximation for d/dx[F(x)] = %.9f and Error Magnitude = %.9f with dx = %.4f \n',i,D,Error_Magnitude(i),dx(i)) %Printing solution to each question part, and setting number of decimal places to display
end

%1.2

plot(dx,Error_Magnitude,'-b'); %Plotting errors against step length
title('Error Magnitude of approximation against $\Delta x$','interpreter','latex'); %Title added to graph
xlim([0 max(dx)]) %Setting x axis limits 
ylim([0 max(Error_Magnitude)]) %Setting y axis limits
xlabel('$\Delta x$','interpreter','latex'); %Labelling x axis
ylabel('Error Magnitude'); %Labelling y axis
legend('Error Magnitude'); %Adding legend to graph

disp('Observe that by increasing the step length by a factor of two we obtain successive magnitude error reductions of factor: ') %Explanatory text
Factor = Error_Magnitude(1:5)./Error_Magnitude(2:6) %Calculating and printing the ratio between respective errors after step length changes
disp('when comparing magnitude errors in (1) & (2), (2) & (3),...,(5) & (6)') %Explanatory text

%From the printed values we see the error magnitude decreases as delta x
%decreases meaning it becomes more accurate, and this is illustrated by the
%plotted graph. Further we can see that the rate of increase with
%delta x is non-linear.
%From the observed values and command window output 'Factor' we see that 
%reducing the step length by a factor of 2 reduces the error magntiude by
%by a factor of 4. With this information we can deduce that the rate of
%increase as step length increases is quadratic.