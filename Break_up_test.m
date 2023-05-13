%testing the SBM from MOCAT

%Russian Asat Test

%Characteristics of Satellite
p1 = [];
p1(1) = 1750; % mass kg (Also seen at 2200 kg)
p1(2) = 2.5; % radius m (from leolabs)
p1(3) = 3;%object class(not a rocketbody)
% p1.r = [0, 0, 0];
% p1.v = [0, 0, 0];

%Characteristics of ASAT (Taken from Raytheon Exoatmospheric Kill Vehicle)
p2 = [];
p2(1)= 63.5029; % mass kg
p2(2) = 1.0033; %radius m (24 in in diam, 55 in length, taking the average for 39.5 in)
p1(3) = 3; %object class(not a rocketbody)
%relative velocity from IAA
dv = 4.6; %km/s (error bounds 3-6 km/s)


%Epoch of event
ep = datetime(2021, 11, 15, 2, 47, 31, 500); % UTC

[fragements1, fragments2] = frag_col_SBM_vec(p1, p2, dv);
tot_KOS1408_frag = length(fragments2(:, 1)) + length(fragements1(:, 1))


%Long March 6A - 350 pieces

%Characteristics of Rocket Body
p1 = [];
p1(1) = 21000; % mass kg 
p1(2) = (3.5+7.3)/2; % radius 5m, height 7.3m
p1(3) = 6;%object class (rocket body)

fragment_LM6A = frag_exp_SBM_vec(p1);
tot_LM6A_frag = length(fragmentsrb(:,1))

%H2-A - 23 pieces 

%Characteristics of Rocket Body
p1 = [];
p1(1) = 20000; % mass kg 
p1(2) = (4+9.2)/2; % radius 4m, height 9.2m
p1(3) = 6;%object class (rocket body)

fragment_H2A = frag_exp_SBM_vec(p1);
tot_H2A_frag = length(fragmentsrb(:,1))

%Orbcomm - 7 Pieces (treating as explosion first, may be collision with small object)

%Characteristics of satellite
p1 = [];
p1(1) = 42; % mass kg 
p1(2) = 3.5; % radius 4m, height 9.2m
p1(3) = 6;%object class (rocket body)

fragment_Orbcomm = frag_exp_SBM_vec(p1);
tot_orbcomm_frag = length(fragmentsrb(:,1))

%Kosmos 2499 - 85 pieces (treating as explosion first, may be collision with small object)

%Characteristics of satellite
p1 = [];
p1(1) = 0; % unknown kg 
p1(2) = 0; % radius unknown
p1(3) = 6;%object class (rocket body)

fragment_Orbcomm = frag_exp_SBM_vec(p1);
tot_orbcomm_frag = length(fragmentsrb(:,1))



