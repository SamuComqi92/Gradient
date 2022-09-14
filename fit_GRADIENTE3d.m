function u = fit_GRADIENTE3d

printf("\nI'm reading the data from the .txt file...\n")

data = load('Data_Astro.txt');          #Reading data
x = data(:,1)-0.8;                      #Eddington ratio corrected
y = data(:,2)-46;                       #Luminosity corrected
z = data(:,3);                          #Luminosity ratio
m = length(y);                          #Number of records
Xe = [ones(m,1),x,y];                   #Add bias unit

#Initial parameters
ITER = 10000;
precision = 1e-06;
alphamin = 0.0001;
alphamax = 10;
n_alpha = 20;
step = (log10(alphamax) - log10(alphamin))/n_alpha;              #Step size for the learning rate

#Loop to find the best learning rate
printf("\nStarting the gradient descent algorithm...")
printf("\nI will use different value of alpha to find the best one with the lowest cost function")
printf("\nThe algorithm will stop when a precision on THETA of 1e-06 is reached!\n")
REStot = [];
alpha = alphamin;
for kk=1:n_alpha,
     theta = [0;0;0];
     RES = [];
     RESULTS2 = [];
     for i=1:ITER,
          ecco = Xe'*(Xe*theta - z)*alpha/m;           #Change
          ecco2 = (sum((Xe*theta-z).^2))/(2*m);        #Cost function
          theta = theta - ecco;                        #Weights update
          RES(i,1) = i;
          RES(i,2) = ecco2;
          RESULTS2(i,:) = ecco(:,1);
          if (abs(RESULTS2(i,1))<=precision) && (abs(RESULTS2(i,2))<=precision) && (abs(RESULTS2(i,3))<=precision),
               break;
          end;
     end;
     REStot(kk,1) = alpha;
     REStot(kk,2) = length(RES);
     REStot(kk,3) = ecco2;
     REStot(kk,4:6) = theta(1:3,1);
     printf("\n[%d/%d] The value of alpha is: %f - The last value of the cost function is: %f (at step %d) \n", kk, n_alpha, alpha, ecco2, length(RES));
     alpha = alpha*(10**(step));
end;

#Gradient descent with the best learning rate
[val, ind] = (min((REStot(:,3))));
theta = [0;0;0];
RES = [];
for i=1:REStot(ind,2),
     ecco = Xe'*(Xe*theta - z)*REStot(ind,1)/m;
     ecco2 = (sum((Xe*theta-z).^2))/(2*m);
     theta = theta - ecco;
     RES(i,1) = i; 
     RES(i,2) = ecco2;
end;

#Comparison with Normal Equation result
printf("\nComparison of the THETA parameters with those calculated with the Normal Equation\n");
C = ((pinv(Xe'*Xe))*Xe'*z)';
last = length(RES);
disp(sprintf("Last iteration number: %.f (precision = %.d)", last, precision))
disp(sprintf("\nPar0_gr = %.4f    Par1_gr = %.4f      Par2_gr = %.4f \nPar0_eq = %.4f    Par1_eq = %.4f     Par2_eq = %.4f \n" , REStot(ind,4), REStot(ind,5), REStot(ind,6), C(1,1), C(1,2), C(1,3)   ))

#Plot of the cost function
printf("\nPlotting the cost function at each iteration\n");
figure(1)
subplot(2,1,1); plot(i=1:last,RES(i,2), 'b', 'Linewidth', 3);
hold on;
grid on;
xlabel("Iterations");
ylabel("J(theta)");
title(sprintf("Best alpha value: %.2f", REStot(ind,1)));
legend(sprintf("Last iteration: %.f", last));

# Procedure to find the parametric equation in 3D space
DATA = [x,y,z]';
MM = mean(DATA,2);
XYZ = DATA - MM;
[U,~,~] = svd(XYZ);
DD = U(:,1);                  #Direction of the parametrix equation
TT = DD'*XYZ;
t1 = min(TT); 
t2 = max(TT);
XYZf = MM + [t1,t2] .* DD;    #Parametric equation
xx = XYZf(1,:);               
yy = XYZf(2,:);
zz = XYZf(3,:);

#Plot data and best-fit function in 3D
printf("\nPlotting the data with the best-fit line\n")
subplot(2,1,2); plot3(x,y,z, 'ro', 'MarkerEdgeColor','k', 'markersize', 3, 'markerfacecolor', 'r')
hold on;
grid on;
plot3(xx,yy,zz, 'b', 'Linewidth', 3)
xlabel("Log Eddington ratio");
ylabel("Observed disk luminosity (10^{46} erg/s)");
zlabel("Observed ratio");
xlim([-2.5 0]);
ylim([-1 1]);
zlim([-2 2]);
view(-10, 10);
legend("Data",sprintf("Best fit   z = %.2f x + %.2fy %.2f", REStot(ind,5), REStot(ind,6), REStot(ind,4) ));

end;
