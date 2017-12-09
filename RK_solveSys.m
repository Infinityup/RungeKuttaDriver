function [W_final,error,file_name] = RK_solveSys(s,t0,final,w0,dt,F,exact)
% Given a a system of IVP's, the desired time step, and the desired number
% of stages, s, this function will solve the system using an s-stage Runge
% Kutta method.
%
% Input:    s = the number of stages for the RK method
%           t0 = the initial time
%           final = the final time
%           w0 = vector of initial conditions for each DE in the system
%           dt = the size of the time step
%           slope = the system of ODE's to be solved
%           exact = (optional) the known exact solution, W(t), used to
%                   calculate error
%
% Output:   W_final = the approximate solution at t = final
%           error = the absolute error at t = final; abs_error = Inf if
%                   the exact solution is not provided
%           file_name = the name of the text file the soln is printed to
%
% Note: all time steps and the approximation at each time step are saved to
%       a .txt file, whose name is in the following format:
%       '(s)stageRK_soln_(step size).txt'
%       e.g. 3-stage RK with step size = 1/4 >> '3stageRK_soln_.2500.txt'

    N = length(w0);
    %initialize output in case of error later on
    W_final = zeros(N,1);
    error = Inf;

    %choose appropriate routine to set coefficients
    [A,b,c] = chooseCoeffs(s);

    %if b is empty, chooseCoeffs ecountered an error and program is ended
    if b == 0
        return
    end

    %open text file for output
    file_name = sprintf('%dstageRK_syssoln_%.4f.txt',s,dt);
    outfile = fopen(file_name,'w');
    line_form = '';
    data_form = '%20.16f\t';
    for i = 1:N+1
        line_form = strcat(line_form,data_form);
    end
    line_form = strcat(line_form,'\n');
    %start with initial condition
    t = t0;
    W = w0;

    %while final time is not yet reached
    while t <= final
        %print t and approx at t
        fprintf(outfile,line_form,t,W);

        %once final is reached, stop the method
        if (final - t) <= dt
            W_final = W;
            break
        end

        %calculate the next approx
        W = RKStep(A,b,c,W,t,dt,F);
        %go forward to the next time step
        t = t + dt;

        %check that the solution is not unbounded, if unbounded, end the
        %program
        if Unbounded(W)
            return
        end

    end

    %close output file
    fclose(outfile);

    %check if exact solution was provided
    if nargin() == 7
        %if yes, calc the abs error at t = final
        w = exact(t);
        error = norm(w - W_final)/norm(w);
    end

end

function [A,b,c] = chooseCoeffs(s)
% Given the number of stages, s, returns the coefficients for an s-stage RK

    switch s
        case 1
            [A,b,c] = oneStage();

        case 2
            [A,b,c] = twoStage();

        case 3
            [A,b,c] = threeStage();

        case 4
            [A,b,c] = fourStage();

        case 6
            [A,b,c] = sixStage();

        otherwise
            fprintf('Error: Invalid input for number of stages.\n')
            fprintf('Program terminating.\n\n')
            A = 0; b = 0; c = 0;
            return

    end
    %check the coefficients satisfy the requirements for RK, if not print
    %error message and return all coefficients as zero
    if ~checkCoeffs(A,b,c)
        fprintf('Error: Coefficients did not satisfy requirements for Runge Kutta.\n')
        fprintf('Program terminating.\n\n')
        A = 0; b = 0; c = 0;
    end
end

function [A,b,c] = oneStage()
% Returns the coefficients for the forward Euler method
    A = 0;
    c = 0;
    b = 1;
end

function [A,b,c] = twoStage()
% Returns the coefficients for the midpoint method
    A = [0   0   ;...
         1/2 0 ];
    b = [0 1]';
    c = [0 1/2]';
end

function [A,b,c] = threeStage()
% Returns the coefficients for a 3-stage RK method
    A = [0  0  0 ;...
        1/2 0  0 ;...
        -1  2  0];
    b = [1/6 2/3 1/6]';
    c = [0 1/2 1]';
end

function [A,b,c] = fourStage()
% Returns the coefficients for a 4-stage RK method
    A = [0  0   0  0 ;...
        1/2 0   0  0 ;...
        0   1/2 0  0 ;...
        0   0   1  0];
    b = [1/6 1/3 1/3 1/6]';
    c = [0 1/2 1/2 1]';
end

function [A,b,c] = sixStage()
% Returns the coefficients for a 6-stage RK method
    A = [0        0          0          0         0      0 ;...
        1/4       0          0          0         0      0 ;...
        3/32      9/32       0          0         0      0 ;...
        1932/2197 -7200/2197 7296/2197  0         0      0 ;...
        439/216   -8         3680/513   -845/4104 0      0 ;...
        -8/27     2          -3544/2565 1859/4104 -11/40 0];
    b = [16/135 0 6656/12825 28561/56430 -9/50 2/55]';
    c = [0 1/4 3/8 12/13 1 1/2]';
end

function [result] = checkCoeffs(A,b,c)
% Checks that all the coefficients satisfy the requirements for RK
    s = length(b);

    %check that the sum of all b's is 1
    if abs(1 - sum(b)) < 1e-15
        b_check = 1;
    else
        b_check = 0;
    end

    %check that c(i) = sum of A(i,j), j = 1:s
    for i = 1:s
        if abs(c(i) - sum(A(i,:))) < 1e-15
            c_check = 1;
            continue
        else
            c_check = 0;
            break
        end
    end

    %if the b_check and the c_check are successful, the coefficients are
    %correct, if either b_check and c_check is unsuccesful the coefficients
    %are incorrect and an error message will be printed by chooseCoeffs()
    if b_check && c_check
        result = 1;
    else
        result = 0;
    end

end

function [W] = RKStep(A,b,c,w0,t,dt,F)
% Returns the approximation at the next time step given the appropriate
% coefficients and the previous time step and the previous approximation.
    s = length(b);
    N = length(w0);
    K = zeros(N,s);


    for i = 1:s
        K(:,i) = dt * F(t + (c(i)*dt) , w0 + (K * A(i,:)'));
    end


    W = w0 + K*b;

end

function [flag] = Unbounded(W)
% Checks that the solution has not become unbounded
    bound = 1e4;

    %if solution is greater than the bound 10000, return 'true' and print
    %an error message
    if (abs(W) >= bound)
        flag = 1;
        fprintf('Solution has become unbounded. \n Program Terminating\n')
    else
        flag = 0;
    end
end
