% Lab 3, #2 : Model of Simple Virus
%set parameters
r = .15;
k = .025;
p = .1;


%set initial conditions
t0 = 0;
a0 = 0;
v0 = .01;
w0 = [a0;v0];
%pick the desired stage for RK
s = 4;

%set final time and time step
final = 96;
dt = .25;

%viral threshold to determine if patient survives or not
viral_thresh = 20;
%color code for patient status
dead = 'r';
alive = 'g';


while p > 1e-6
    %the system to be solved
    F = @(t,W) [(k*W(2)) ;( r*W(2) - p*W(2)*W(1) )];
    
    %solve the system
    [W,e,data_file] = RK_solveSys(s,t0,final,w0,dt,F);
    
    %read the solution from file
    infile = fopen(data_file,'r');
    data = fscanf(infile,'%f%f%f',[3,Inf])';
    fclose(infile);
    
    %part (a): plot both densities for p = .1
    if p == .1
        %plot both the virus density and the antibody density
        
        subplot(1,2,1)
        plot(data(:,1),data(:,2),'g')
        hold on
        plot(data(:,1),data(:,3),'b')
        title('Virus and Antibody Density over 4 Days')
        xlabel('time (minutes)')
        ylabel('Virus/Antibody density')
        legend('Antibody','Virus')
        hold off
    else
        hold on
    end
    
    %part (b) plot jut the viral denity as p is decreased
    subplot(1,2,2)
    %check if patient has died and plot accordingly
    if sum(data(:,3) >= viral_thresh)
        %if patient has died, plot and then terminate the model
        plot(data(:,1),data(:,3),dead)
        fprintf('Patient has died. Viral threshold of %f reached. \n',viral_thresh)
        fprintf('p(t) = %f at event of patient death \n',p)
        break
    else
        plot(data(:,1),data(:,3),alive)
    end
    
    %decrease p by .01
    p = p - .01;
    

end

%label the second plot
xlabel('time (minutes)')
ylabel('Virus density')
title('Virus density with different p(t) values')
hold off

