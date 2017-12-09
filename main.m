%% ISC 4232 : Lab Assignment #3
% Lydia LaSeur
% October 7, 2014

%% Part 1: Modify Existing RK Library
% Modify the previously completed explicit Runge Kutta library to handle a
% system of firt order IVPs.  Run and test the convergence rates for stages
% 1 through 4 using the driver and the following system of three first
% order IVPs.
%%
% w1'(t) = 2*w2(t) - 4t
%%
% w2'(t) = -w1(t) + w3(t) - e^t +2
%%
% w3'(t) = w1(t) - 2*w2(t) + w3(t) + 4t
%%
% w1(0) = -1 ; w2(0) = 0 ; w3(0) = 2 ; 0 < t < 10

RK_sysDriver

%% Part 2: Model of Simple Viral Infection
% (a). Model a simple viral infection using the following system and given
% parameters.  Plot the density of the virus, v(t), and the density of
% antibodies, a(t) over the given domain, in time t.
%
%%
% a'(t) = k*v(t)
%%
% v'(t) = r*v(t) - p*v(t)a(t)
%%
% a(0) = 0 ; v(0) = .01 ; 0 < t < 96 (hours)
%%
% k = .025 ; r = .15 ; p = .1
%%
% Where k is the proportionality constant (number of antibodies produced is
% proportional to the number of viruses present).  r is the rate at which
% new viruses occur.  p is the probability that at any gven instance a
% particular antibody will encounted a particular virus and successfully
% destroy it.  See both solutions on the first plot below.

%%
% (b). Decrease the value of p by .01 and run the model for each new value
% of p.  PLot the solution of the viral denity for each run to determine
% when the patient will die, i.e. when the viral density exceeds 20. Keep
% all other parameters and initial conditions the same as before.  See
% results in the second plot below.

VirusModel

%%
% The color of the curve correponds to the status of the patient.  Green
% means the patient survived and red means the patient died.

%% Part 3: Model Viral Infection with Single Mutation
% Modify the previous model to account for mutation of the virus.  To keep
% it simple, the virus will only mutate once.  Since the virus will mutate
% different antibodies will be produced to fight the new mutated virus.
% This means two more DE mut be added to the system, one for the viral
% density of the mutated strain v2'(t), and one for the density of the
% antibodies that will fight the mutated virus a2'(t).  Plus, a DE the
% describes the general immune response g'(t), i.e. production of
% antibodies to fight any virus, not a specific one.  The 'strength' of the
% the general response is denoted by a constant, s.  The virus is also able
% fight the antibodies.  The strength of the virus is denoted by the
% constant u.  The model is described by the five equation system below.
%%
% v1'(t) = v1(t) * (r - s*g(t) - p*a1(t)
%%
% a1'(t) = k*v1(t) - u*a1(t) * (v1(t) + v2(t))
%%
% v2'(t) = v2(t) * (r - s*g(t) - p*a2(t)
%%
% a2'(t) = k*v2(t) - u*a2(t) * (v1(t) + v2(t))
%%
% g'(t) = (k - u*g(t)) * (v1(t) + v2(t))
%%
% a1(0) = a2(0) = 0 ; v1(0) = .01 ; v2(0) = 0 ; g(0) = 0 ; 0 < t < 100

VirusModelOneMutation
%%
% The plot above shows the total viral density and the total antibody
% density.  The first big peak in the viral density is the mostly the
% original, non-mutated virus.  At 10 hours, the virus mutates and second
% virus is introduced into the system.  As the a1 antibodies fight the
% original virus the viral density begins the decrease but now the density
% of the mutated strain increases rapidly since the a2 antibodies haven't
% had time to replicate.  This is the second peak in the viral density.
% The viral density will continue to oscillate in this fashion, while
% continuously decreasing. While the antibody sum never decreases the plot
% below shows that the antibody denisty for the specific anitboies, a1 and
% a1, do oscillate (due to the addition of the virus fighting back)but in
% a manner that they're behavoir mirrors each other, enuring that their sum
% doesn't decrese.  In our model, the virus cannot fight the general
% antibodies, g, so that denisty never decreases.  In this model using
% these parameters, the patient survives since the viral density does not
% exceed 20.

%plot the inidiviudal antibody densities
plot(soln(:,1),soln(:,3),'b') %a1(t) density
hold on
plot(soln(:,1),soln(:,5),'r') %a2(t) density
plot(soln(:,1),soln(:,6),'g') %g(t) desnity
xlabel('Time (hours)')
ylabel('Antibody Density')
title('Density of a1, a2 and g Antibodies')
legend('a1','a2','g')
