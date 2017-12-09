%test stages 1 through 4 in RK_solveSys
stages = [1 2 3 4];
dt = .1;
step_sizes = dt * [1 1/2 1/4 1/8 1/16]';

%IVP system and exact solution to system
t0 = 0;
final = 10;
w0 = [-1 0 2]';

F = @(t,W) [(2*W(2) - 4*t)               ;...
            (-W(1) + W(3) - exp(t) + 2)  ;...
            (W(1) - 2*W(2) + W(3) + 4*t)];
        
exact = @(t,W)[(-cos(2*t))         ;...
               (sin(2*t) + 2*t)    ;...
               (cos(2*t) + exp(t))];
           
N = length(w0);
%for each method
for i = 1:length(stages)
    %print table header
    fprintf('%d-STAGE RUNGE KUTTA METHOD FOR %d EQN SYSTEM     \n',stages(i),N)
    fprintf('Step Size    Abs Error at T      Convergence Rate \n')
    fprintf('---------    --------------      ---------------- \n')
    errs = zeros(length(step_sizes),1);
    
    %for each step size
    for j = 1:length(step_sizes)
        %solve the IVP using the current method and current step size
        [soln,errs(j)] = RK_solveSys(stages(i),t0,final,w0,step_sizes(j),F,exact);
        
        %print the current step size, current error and convergence rate to
        %table
        if j == 1
            fprintf('%-13.6f%-20.16f n/a \n',step_sizes(j),errs(j))
        else
            r = log(errs(j)/errs(j-1)) / log(step_sizes(j)/step_sizes(j-1));
            fprintf('%-13.6f%-20.16f%-13.12f \n',step_sizes(j),errs(j),r);
        end
    end
    fprintf('------------------END OF TABLE------------------\n\n\n')
    
end

