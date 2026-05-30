%__________________________________________________________________
%  Heavy Ball Optimizer (HBO)
%  Developed in MATLAB R2024a

%  programmer: Mohammed Jameel  
%  E-Mail: moh.jameel@su.edu.ye, mohjameel555@gmail.com
 
%  Paper:
%  Mohammed Jameel, Ahmed R. El-Saeed, Anis Elgabli, Ibrahim Al-Shourbaji, Pramod H. Kachare, Mohamed Abouhawwash 
%  Heavy Ball Optimizer: A Momentum-Driven Metaheuristic with Superior 
% Scalability and Real-World Applications          
% _________________________________________________________________________%
% ---------------------------------------------------------------------------------------------------------------------- %

clear;
close all;
clc;

Population=50;      % Number of Population

Func_name='F1';     % Name of the Test Function

MaxIt=800;          % Maximum number of iterations

% Load details of the selected benchmark function
                         [LB,UB,dim,fobj]=Get_Functions_details(Func_name);

  [Best_Score,Best_Pos,Convergence_curve] = HBO(Population,MaxIt,LB,UB,dim,fobj);
%% Plots
figure('Position',[500 500 660 290])

% Draw Search Space
subplot(1,2,1);
func_plot(Func_name);
title('Parameter Space')
xlabel('X_1');
ylabel('X_2');
zlabel([Func_name,'( X_1 , X_2 )'])

% Draw Objective Space
subplot(1,2,2);
semilogy(Convergence_curve,'Color','[1 0 1]','LineWidth', 1.5)
title('Objective Space')
xlabel('Iteration');
ylabel('Best Fitness');
axis tight
grid on
box on
legend('HBO')

display(['Best solution obtained by HBO is : ', num2str(Best_Pos)]);
display(['Best fitness found by HBO is     : ', num2str(Best_Score)]);





