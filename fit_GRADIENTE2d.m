function u = fit_GRADIENTE2d

printf("\nI'm reading the data from the .txt file...\n")

data = load('Data_Astro.txt');       #Load data
x = data(:,1)-0.8;                   #Eddington Ratio corrected
y = data(:,3);                       #Luminosity ratio
m = length(y);                       #Number of records
X1 = [ones(m,1),x];                  #Add bias unit

#Initial parameters
ITER = 10000;
precision = 1e-06;
alphamin = 0.0001;
alphamax = 10;
n_alpha = 20;
step = (log10(alphamax) - log10(alphamin))/n_alpha;                    #Step size for the learning rate

#Loop to find the best learning rate
printf("\nStarting the gradient descent algorithm...")
printf("\nI will use different value of alpha to find the best one with the lowest cost function")
printf("\nThe algorithm will stop when a precision on THETA of 1e-06 is reached!\n")
REStot = [];
alpha = alphamin;
for kk=1:n_alpha,
     theta = [0;0];
     RES = [];
     RESULTS2 = [];
     for i=1:ITER,
          ecco = X1'*(X1*theta - y)*alpha/m;                   #Change
          ecco2 = (sum((X1*theta-y).^2))/(2*m);                #Cost function
          theta = theta - ecco;                                #Weights update
          RES(i,1) = i; 
          RES(i,2) = ecco2;
          RESULTS2(i,:) = ecco(:,1);
          if (abs(RESULTS2(i,1))<=precision) && (abs(RESULTS2(i,2))<=precision),
               break;
          end;
     end;
     REStot(kk,1) = alpha;
     REStot(kk,2) = length(RES);
     REStot(kk,3) = ecco2;
     REStot(kk,4:5) = theta(1:2,1);
     printf("\n[%d/%d] The value of alpha is: %f - The last value of the cost function is: %f (at step %d) \n", kk, n_alpha, alpha, ecco2, length(RES));
     alpha = alpha*(10**(step));

end;

#Gradient descent with best learning rate
[val, ind] = (min((REStot(:,3))));
theta = [0;0];
RES = [];
for i=1:REStot(ind,2),
     ecco = X1'*(X1*theta - y)*REStot(ind,1)/m;
     ecco2 = (sum((X1*theta-y).^2))/(2*m);
     theta = theta - ecco;
     RES(i,1) = i; 
     RES(i,2) = ecco2;
end;

#Comparison with Normal equation results
printf("\nComparison of the THETA parameters with those calculated with the Normal Equation\n");
C = ((pinv(X1'*X1))*X1'*y)';
last = length(RES);
disp(sprintf("Last iteration number: %.f (precision = %.d)", REStot(ind,2), precision))
disp(sprintf("\nPar0_gr = %.4f    Par1_gr = %.4f \nPar0_eq = %.4f    Par1_eq = %.4f\n" , REStot(ind,4), REStot(ind,5), C(1,1), C(1,2)  ))


#Plot cost function
printf("\nPlotting the cost function at each iteration\n");
figure(1)
subplot(2,1,1); plot(i=1:last,RES(i,2), 'b', 'Linewidth', 3);
hold on;
grid on;
xlabel("Iterations");
ylabel("J(theta)");
title(sprintf("Alpha value: %.2f", REStot(ind,1)));
legend(sprintf("Last iteration: %.f", last));
t1 = -2.5:((max(x)-min(x))/100.):0;
F1 = REStot(ind,4) + REStot(ind,5)*t1;

#Plot of data
printf("\nPlotting the data with the best-fit line\n")
subplot(2,1,2); plot(x,y, 'ro', 'MarkerEdgeColor','k', 'markersize', 3, 'markerfacecolor', 'r')
hold on;
grid on;
plot(t1,F1, 'b', 'Linewidth', 3)
xlabel("Log Eddington ratio");
ylabel("Log Observed Ratio");
xlim([-2.5 0]);
ylim([-2 2]);
legend("Data",sprintf("Best fit   y = %.2f x %.2f", REStot(ind,5), REStot(ind,4)));


end;

