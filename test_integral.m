clear all %Clearing all prior variables
close all %Closing all prior graph
format short %Formatting output where not specified

%2.1 & 2.2
x=zeros(5,1); %Creating a vector of length 5 to loop x values into
for i=1:length(x) %Looping over length of x to obtain array of 5 points
    x(i) = 0.5 - 0.5*cos(((2*i-1)*pi)/10); %Looping corresponding Chebychev points into vector using equation (4)
end

weight = zeros(1,length(x)); %Preallocating array to append weights into, of same length as loop, for speed
for i=1:length(x) %Looping 5 times for each L_i(x)
    c=zeros(5,1); %Creating vector of zeros of length 5
    c(i)=1; %Entering 1 to row corresponding to L_i(x_i) from (5)
    coefficients = polyfit(x,c,4); %Finding coefficients to system at L_i
    %I didn't know if the question implied it wanted these displayed or not.
    fprintf('The corresponding x^4,...,x^0 coefficients to L%u(x) are :',i) %Producing text to explain the solution to user
    disp(coefficients) %Displaying solved coefficient system for that L_i
    pint = polyint(coefficients); %Integrating generated polynomial 
    weight(i) = polyval(pint,1)-polyval(pint,0); %Calculating weights for each point and appending to weight array - concept taken from lab 7
end

disp('The weights W_1,...,W_5 are :');%Explanation of following output
disp(weight); %Printing weights corresponding to explanation above

%2.3
k=0; %Stating initial degree of polynomial
syms x_symbol; %Defining x for symbolic integration - taken from: https://uk.mathworks.com/help/symbolic/integration.html?s_tid=gn_loc_drop
f = x_symbol; %Defining function to be integrated to some powers of k
integral = double(int(f^k,0,1)); %Integrating f=x^0 over [0,1]
approximation = sum(weight.*(x'.^k)); %Calculating approximate integral for x^0, together lines 30,31 test whether it is already degree 0 precision

initial_error = abs(integral-approximation); %Calculting original error
error = [initial_error]; %Creating array to loop magnitude error values into, starting with initial value

%Note to examiner - on my laptop (MATLAB 2017a) 1e-15 - MATLAB standard dp - works as the first approximation above is stored as 1.000000000000000, while on the library computer it
%stores the first approximation as 1.000000000000004 and so skips the loop - I believe this is just a bit/version issue https://uk.mathworks.com/matlabcentral/answers/37979-different-results-when-running-the-same-matlab-program-on-different-version-of-matlab
%hence 1e-14 is used as a precaution.
while abs(integral - approximation)<1e-14 %Running loop until the error is no longer attributed to round-off (round-off error happens at 15(library pc)/16(laptop) decimal places)
    k = k+1; %Increasing degree for test
    integral = double(int(f^k,0,1)); %Integrating f=x^k between [0,1]
    approximation = sum(weight.*(x'.^k)); %Calculating approximate integral for x^k.
    error(k+1)=abs(integral-approximation); %Appending errors for each corresponding k, starting at second position
end

if k <= 0 %Testing if the model never reflects the real integral due to round off error
    fprintf('Degree of precision never sufficient for this integral\n') %Printing error
else %Otherwise output degree of precision
    fprintf('The quadrature rule has degree of precision %u. \n',k-1) %Printing out degree of precision, that is, the last value k for which the error was less than round-off.
    x_plot = 0:k; %Creating array of x-axis points to plot corresponding to order of polynomials
    plot(x_plot,error,'-r'); %Plotting the error magnitude against order of the polynomial
    legend('Error Magnitude') %Creating legend
    xlim([0 max(x_plot)]); %Setting x axis limits 
    ylim([0 max(error)]); %Setting y axis limits
    xlabel('Order of Polynomial'); %Labelling x axis
    ylabel('Magnitude Error'); %Labelling y axis
end

%This output demonstrates that this weight method calculates integrals up
%to degree 5 accurately beyond round-off error. This is illustrated by the
%figure which demonstrates a jump in error at n=6 to more than 4.5*10^(-5)
%which is a change in magnitude error greater than could be considered simply round-off error,
%with a magnitude error of more than (10)^-5. The theoretical expected
%value for a Gaussian Quadrature is accurate for polynomials of degree at most 2n-1 for n points, 
%so in our case 5<2(5)-1=9. Hence, the result is not unexpected but less
%accurate implying a less suitable choice of points and weights.