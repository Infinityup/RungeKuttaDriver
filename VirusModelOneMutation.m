% Lab 3, #3 : Model Virus with One Mutation
%set parameters
k = .1;
r = .5;
p = .5;
s = .25;
u = .1;

%set inital conditions
a1_0 = 0;
a2_0 = 0;
v1_0 = .01;
v2_0 = 0;
g0 = 0;
W0 = [v1_0, a1_0, v2_0, a2_0, g0]';
t0 = 0;

%set the time at which virus will mutate and time step
mutate = 10;
dt = .25;

%set the desired stage for RK
stage = 4;


%the system to be solved
F = @(t,W) [W(1)*(r - s*W(5) - p*W(2))     ; ...
            k*W(1) - (u*W(2) * (W(1)+W(3))); ...
            W(3)*(r - s*W(5) - p*W(4))     ; ...
            k*W(3) - (u*W(4) * (W(1)+W(3))); ...
            (k - u*W(5)) * (W(1)+W(3)) ]   ;


%solve the system until t = 10 hours      
[W_10,e,data_file] = RK_solveSys(stage,t0,mutate,W0,dt,F);
%read in the solution from file
infile = fopen(data_file,'r');
soln = fscanf(infile,'%f%f%f%f%f%f',[6,Inf])';
fclose(infile);

%mutate the virus
W_10(3) = .0001;
%set mutate time as the initial time and set the final time
t0 = mutate;
final = 100;
%solve the system again
[final,e,data_file] = RK_solveSys(stage,t0,final,W_10,dt,F);
%read in the solution and add to the previous solution
infile = fopen(data_file,'r');
soln = [soln;fscanf(infile,'%f%f%f%f%f%f',[6,Inf])'];

%calc the viral and antibody sums
virus_sum = soln(:,2) + soln(:,4);
antib_sum = soln(:,3) + soln(:,5) + soln(:,6);

%plot both sums
subplot(1,1,1)
plot(soln(:,1),virus_sum,'b');
hold on
plot(soln(:,1),antib_sum,'g');
xlabel('Time (hours)')
ylabel('Virus/Antibody Density')
title('Total Viral and Antibody Density (One Mutation)')
legend('Viral Sum','Antibody Sum')
hold off



