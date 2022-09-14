% CODE BY: Camino Rodríguez Pérez-Carral, Rodrigo de Lera de las Heras and Jorge Lázaro Ruiz
% USAGE OF THIS CODE IS STRICTLY FOR REFERENCE ONLY, DO NOT COPY

% Rodrigo de Lera de las Heras 100452323
% Jorge Lázaro Ruiz 100452172
% Camino Rodríguez Pérez-Carral 100445091

clear;
clc;

%contour plot of the Himmelblau’s function
x = linspace(-6,6);
y = linspace(-6,6);%in the square[-6,6]
[X,Y] = meshgrid(x,y);
Z = (X.^2 + Y -11).^2 + (X+Y.^2 -7).^2;
contour(X,Y, Z,100)
title("Newton's method in the Himmelblau’s function");
hold on

%grid search
fprintf('Wait a second, we are working on the grid search...\n');

for i = -6:0.12:6
    for j=-6:0.12:6
        pointforGS = [i,j]; %get a point i, j
        xGridSearch = functionGridSearch(pointforGS); %the functionGridSearch takes the x of the point it converges to
        xGridSearch = round(xGridSearch, 6); %we round it to 6 decimals so that there is no problem comparing it with the extrema
        switch xGridSearch
            %minima:
            case 3
                hold on
                plot(i,j, '--.r')
            case -2.805118
                hold on
                plot(i,j, '--.b')
            case -3.779310
                hold on
                plot(i,j, '--.g')
            case 3.584428
                hold on
                plot(i,j, '--.y')
            %maxima
            case -0.270845
                hold on
                plot(i,j, '--.m')
            %saddle points
            case {-3.073026, 3.385154, -0.127961, 0.086678}
                hold on
                plot(i,j, '--.c')
        end
    end
end
hold off

%points
x1 = [7,7];
x2 = [27.18,0];
x3 = [0.628,318];
%we ask the user to select a point
point = input('Select a point: \n1.(-6,-6), \n2.(0.04, 3.09) \n3.(3.01, 3.03)\n4. Exit\n');
while point ~= 4 %goes on and on until exit is pressed
    %and then work on the three cases
    switch point
        case 1
            v1 = funcionNewton(x1);
            fprintf('(a) How many iterations does it take for the method to stop/converge?\n');
            whatIsIt(v1(1),v1(2),v1(3));
            fprintf('(b) What is the condition number of the Jacobian matrix at each step?\n');
            printCd(v1);
        case 2
            v2 = funcionNewton(x2);
            fprintf('(a) How many iterations does it take for the method to stop/converge?\n');
            whatIsIt(v2(1),v2(2),v2(3));
            fprintf('(b) What is the condition number of the Jacobian matrix at each step?\n');
            printCd(v2);
        case 3
            v3 = funcionNewton(x3);
            fprintf('(a) How many iterations does it take for the method to stop/converge?\n');
            whatIsIt(v3(1),v3(2),v3(3));
            fprintf('(b) What is the condition number of the Jacobian matrix at each step?\n');
            printCd(v3);
    end
    point = input('Select a point: \n1.(-6,-6), \n2.(0.04, 3.09) \n3.(3.01, 3.03)\n4. See answers to questions and exit\n');
end

%we answer the rest of the questions
fprintf('From your observations about these 3 points, is the process well-conditioned as a whole? \nWe have seen that the condition numbers for these points were always\nvery close to 1. Based on these three, it is well-conditioned.\n');

%setting up some tables that will help illustrate our conclusions
        StartingPoint = ["(-7,7)";"(0,-1000)";"(-2000,-9000)";"(-1618,0)";"(27.18,0)";"(0.628,318)";"(7,7)"]
        ConvergesTo = ["(-2.805118,3.131313)";"(3.584428,-1.848127)";"(-3.779310,-3.283186)";"(-2.805118,3.131313)";"(3.584428,-1.848127)";"(-2.805118,3.131313)";"(3,2)"];
        IterationsNeeded = ["8";"21";"25";"21";"11";"17";"8"];
        Condition = ["Well-conditioned";"Ill-conditioned";"Ill-conditioned";"Ill-conditioned";"Ill-conditioned";"Ill-conditioned";"Ill-conditioned"]
        
        tOutside = table(StartingPoint,ConvergesTo,IterationsNeeded,Condition);
fprintf('(1) What happens if the point’s orbit does, at some step, leave the square [−6, 6] × [−6, 6]? \n       To test this, we placed the initial point outside of this range.\nNo matter how far we placed it, the method always converged to a minimum. The only\ndifference was that it took more iterations, which was to be expected as some of these initial\napproximations were very off.\n');
disp(tOutside)
fprintf('(2) What happens near the crest of the hill in the middle of the square? (Observe the plotted function) \n       Looking at the results of the grid search, we see that the points closer to the hill converge to\nthe maximum point of the function. This can be expected since we are using Newton’s method\nwith an approximation very close to this critical point. If we move away from the crest of the\nhill, the points start to converge to the saddle, as can be seen in the figure.\n');
fprintf('(3) What criteria did you choose for stopping? \nWe decided to use the distance to a given value as our stopping criterion, in this case, we\n established that the difference between the value of the function and 0 is less than 10^-6. We\nchose this number because the minima in the Himmelblau function are usually represented\nwith six decimal figures, so we thought matching this precision would be correct.\nAt first, we also thought about adding a stopping criterion based on a maximum\nnumber of iterations, but after trying our algorithm with a big number of points, it rarely\nsurpassed 20 iterations, so we concluded that establishing a maximum number of iterations\nwas not necessary.');

%functions

function[vfinal] = funcionNewton(x) %recieves a vector containing the x and y of a point, and applies newton's method to it.
%It returns a vector with the point it converges to, the number of iterations, and the conditions numbers
fprintf("Applying Newton-Raphson to (%.6f,%.6f)\n", x(1), x(2));
%we made the calculations by hand because this way the program is a lot faster
f = [4*x(1)*((x(1)^2)+x(2)-11)+2*(x(1)+(x(2)^2)-7),2*((x(1)^2)+x(2)-11)+4*x(2)*(x(1)+(x(2)^2)-7)]; %Gradient of the Himmelblau function
counter = 0; %This will count the iterations needed
while (abs(f(1))+abs(f(2))>10^(-6)) %Our tolerance will be 10^-6
    f = [4*x(1)*((x(1)^2)+x(2)-11)+2*(x(1)+(x(2)^2)-7),2*((x(1)^2)+x(2)-11)+4*x(2)*(x(1)+(x(2)^2)-7)];
    J = [12*(x(1)^2)+4*x(2)-42,4*x(1)+4*x(2);4*x(1)+4*x(2),12*(x(2)^2)+4*x(1)-26]; %Jacobian matrix of f
    M = inv(J);
    prevF = f; %This is f(x^k-1)
    prev = x; %This is x^k-1
    x = transpose(prev) - M*(transpose(prevF)); % x^k = (x^k-1) -  J^-1(x^k-1) * f(x^k-1)
    x = transpose(x); %In the previous step we needed to transpose the vectors in order to perform correct matrix multiplication, and this line of code reverts it to horizontal form
    counter = counter +1; %update the counter
    
    %now we calculate the condition number of the jacobian matrix
    eigJ = eig(J); %calcualte the eigenvalues
    eigJ = abs(eigJ); %vector of eigenvalues in absolute value
    k = (max(eigJ(1),eigJ(2)))/(min(eigJ(1),eigJ(2))); %calculate the condition number k
    vk(counter) = k; %store that value in a vector
end
%we create and fill the vector it returns
vfinal(1) = counter; %first the iterations
vfinal(2) = x(1); %then the x of the point it converges to
vfinal(3) = x(2); %and the y
for i=1:counter
    vfinal(i+3)=vk(i); %and then fill it with the condition number at each step
end
end

function whatIsIt(i,x,y) %it receives a point(the point where it converged) and the number of iterations and datermines and print what that point is
    z = round((x^2+y-11)^2 + (x+y^2-7)^2,6); % Due to round-off errors, z is not exactly 0, we fix this by rounding it to 6 decimal places (the same precision has been used everywhere else in the exercise);
    if z == 0 %the minimum points of the Himmelblau function have z = 0
        fprintf("The method converges to the global minimum (%.6f,%.6f) in %d iterations.\n", x, y, i);
    elseif round(z,4) == 181.6165 %this is the only local maximum of the function. It is rounded to 4 decimal places because that is the precision where we stop having roundoff errors
        fprintf("The method converges to the local maximum (%.6f,%.6f) in %d iterations.\n", x, y, i);
    else %it is a saddle
        fprintf("The method converges to (%.6f,%.6f) in %d iterations, but the point is not an extremum, it is a saddle.\n", x, y, i);
    end
end

function[result] = functionGridSearch(x) %recieves a vector containing the x and y of a point, and applies newton's method to it.
%It returns the x of the point it converges to. It is a simpler version of
%the functionNewton
f = [4*x(1)*((x(1)^2)+x(2)-11)+2*(x(1)+(x(2)^2)-7),2*((x(1)^2)+x(2)-11)+4*x(2)*(x(1)+(x(2)^2)-7)]; %Gradient of the Himmelblau function
while (abs(f(1))+abs(f(2))>10^(-6)) %Our tolerance will be 10^-6
    f = [4*x(1)*((x(1)^2)+x(2)-11)+2*(x(1)+(x(2)^2)-7),2*((x(1)^2)+x(2)-11)+4*x(2)*(x(1)+(x(2)^2)-7)];
    J = [12*(x(1)^2)+4*x(2)-42,4*x(1)+4*x(2);4*x(1)+4*x(2),12*(x(2)^2)+4*x(1)-26];
    M = inv(J);
    prevF = f; %This is f(x^k-1)
    prev = x; %This is x^k-1
    x = transpose(prev) - M*(transpose(prevF)); % x^k = (x^k-1) -  J^-1(x^k-1) * f(x^k-1)
    x = transpose(x); %In the previous step we needed to transpose or vectors in order to perform correct matrix multiplication, and this line of code reverts it to horizontal form
end
result = x(1); %we return the x of the point it converges to
end

function printCd(vect) %it prints the condition numbers that were stored in a vector
    len = length(vect);
    disp('Condition number at each iteration: ');
    for i = 4:len
        fprintf('%.6f, ', vect(i));
    end
    if max(vect(4:len)) > 1.5 %1.5 chosen arbitrarily, it seems like a good place to draw the line between a good and a bad condition number for this exercise
        fprintf('\nThis method is ill-conditioned for this point.\nThe worst condition number is %.6f, and the best one is %.6f', max(vect(4:len)), min(vect(4:len)));
    else
        fprintf('\nThis method is well-conditioned for this point.\nThe worst condition number is %.6f, and the best one is %.6f', max(vect(4:len)), min(vect(4:len)));
    end
    fprintf('\n');
end

