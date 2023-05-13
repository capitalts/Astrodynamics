% 16.346 Astrodynamics Final Project
% Jackie and Tory Smith
% 16 May 2023
clear; close all; clc

%% COSMOS 1408 ASAT test

% Characteristics of Satellite
p1 = [];
p1(1) = 1750;       % mass kg (Also seen at 2200 kg)
p1(2) = 2.5;        % radius m (from leolabs)
p1(3) = 3;          % object class(not a rocketbody)

% Characteristics of ASAT (Taken from Raytheon Exoatmospheric Kill Vehicle)
p2 = [];
p2(1) = 63.5;       % mass kg
p2(2) = 1.00;       % radius m (24 in in diam, 55 in length, taking the average for 39.5 in)
p1(3) = 3;          % object class(not a rocketbody)
% relative velocity from IAA
dv = 4.6;           % km/s (error bounds 3-6 km/s)

% Epoch of event
ep = datetime(2021, 11, 15, 2, 47, 31, 500); % UTC

LB = 0.10;
[fragments1, fragments2] = frag_col_SBM_vec(p1, p2, dv, LB);
tot_KOS1408_frag_10cm = length(fragments1(:, 1));
tot_KOS1408_frag_1m = length(fragments2(:, 1));

% LB = 0.02;
% [fragments1, fragments2] = frag_col_SBM_vec(p1, p2, dv, LB);
% tot_KOS1408_frag_2cm =  length(fragments1(:, 1));

% KOS1408_2to10 = (tot_KOS1408_frag_2cm-tot_KOS1408_frag_10cm);

disp("COSMOS 1408 ASAT Test")
% disp("Total pieces 2-10cm: " + num2str(KOS1408_2to10))
disp("Total pieces 10cm-1m: " + num2str(tot_KOS1408_frag_10cm))
disp("Total pieces >1m: " + num2str(tot_KOS1408_frag_1m))


%% Long March 6A - 350 pieces

% Characteristics of Rocket Body
p1 = [];
p1(1) = 952;        % mass kg 
p1(2) = 3.0;        % radius m
p1(3) = 6;          % object class (rocket body)

LB = 0.10;

fragment_LM6A = frag_exp_SBM_vec(p1,LB);
tot_LM6A_frag = length(fragment_LM6A(:,1));
tot_LM6A_frag_1m = length(find(fragment_LM6A(:,1)>1));

disp("Chang Zheng 6A")
disp("Total pieces 10cm-1m: " + num2str(tot_LM6A_frag))
disp("Total pieces >1m: " + num2str(tot_LM6A_frag_1m))

%% H2-A - 23 pieces 

%Characteristics of Rocket Body
p1 = [];
p1(1) = 1400;       % mass kg 
p1(2) = 8;          % radius m
p1(3) = 6;          % object class (rocket body)

LB = 0.10;

fragment_H2A = frag_exp_SBM_vec(p1,LB);
tot_H2A_frag = length(fragment_H2A(:,1));

tot_H2A_frag_1m = length(find(fragment_H2A(:,1)>1));

disp("HII-A")
disp("Total pieces 10cm-1m: " + num2str(tot_H2A_frag))
disp("Total pieces >1m: " + num2str(tot_H2A_frag_1m))

%% Kosmos 2499 - 85 pieces (treating as explosion first, may be collision with small object)

%Characteristics of satellite
p1 = [];
p1(1) = 50;         % unknown kg 
p1(2) = 3.5;        % radius unknown
p1(3) = 3;          % object class (not rocket body)

% Explosion
LB = 0.10;

fragment_cosmos = frag_exp_SBM_vec(p1, LB);
tot_cosmos_frag = length(fragment_cosmos(:,1));
tot_cosmos_frag_1m = length(find(fragment_cosmos(:,1)>1));


disp("COSMOS 2499 Explosion")
disp("Total pieces 10cm-1m: " + num2str(tot_cosmos_frag))
disp("Total pieces >1m: " + num2str(tot_cosmos_frag_1m))

% Collision
% Characteristics of Secondary 
p2 = [];
p2(1) = 1.0;       % mass kg
p2(2) = 0.10;       % radius m 
p1(3) = 3;          % object class(not a rocketbody)
% relative velocity from IAA
dv = 5;             % km/s 

LB = 0.10;

[fragments1, fragments2] = frag_col_SBM_vec(p1, p2, dv, LB);
tot_KOS2499_frag_10cm = length(fragments1(:, 1));
tot_KOS2499_frag_1m = length(fragments2(:, 1));

% LB = 0.02;
% [fragments1, fragments2] = frag_col_SBM_vec(p1, p2, dv, LB);
% tot_KOS2499_frag_2cm =  length(fragments1(:, 1));

% KOS2499_2to10 = (tot_KOS2499_frag_2cm-tot_KOS2499_frag_10cm);
disp("COSMOS 2499 Collision")

% disp("Total pieces 2-10cm: " + num2str(KOS2499_2to10))
disp("Total pieces 10cm-1m: " + num2str(tot_KOS2499_frag_10cm))
disp("Total pieces >1m: " + num2str(tot_KOS2499_frag_1m))


%% Orbcomm - 7 Pieces (treating as explosion first, may be collision with small object)

%Characteristics of satellite
p1 = [];
p1(1) = 42;         % mass kg 
p1(2) = 3.5;        % radius 4m, height 9.2m
p1(3) = 3;          % object class (not rocket body)

LB = 0.10;

fragment_Orbcomm = frag_exp_SBM_vec(p1, LB);
tot_orbcomm_frag = length(fragment_Orbcomm(:,1));
tot_Orbcomm_frag_1m = length(find(fragment_Orbcomm(:,1)>1));

disp("Orbcomm Explosion")
disp("Total pieces 10cm-1m: " + num2str(tot_orbcomm_frag))
disp("Total pieces >1m: " + num2str(tot_Orbcomm_frag_1m))

% Collision
% Characteristics of Secondary 
p2 = [];
p2(1) = 1.0;       % mass kg
p2(2) = 0.10;       % radius m 
p1(3) = 3;          % object class(not a rocketbody)
% relative velocity from IAA
dv = 5;             % km/s 

LB = 0.10;

[fragments1, fragments2] = frag_col_SBM_vec(p1, p2, dv, LB);
tot_Orbcomm_frag_10cm = length(fragments1(:, 1));
tot_Orbcomm_frag_1m = length(fragments2(:, 1));

% LB = 0.02;
% [fragments1, fragments2] = frag_col_SBM_vec(p1, p2, dv, LB);
% tot_Orbcomm_frag_2cm =  length(fragments1(:, 1));

% Orbcomm_2to10 = (tot_Orbcomm_frag_2cm-tot_Orbcomm_frag_10cm);

% disp("Total pieces 2-10cm: " + num2str(Orbcomm_2to10))
disp("Orbcomm Collision")

disp("Total pieces 10cm-1m: " + num2str(tot_Orbcomm_frag_10cm))
disp("Total pieces >1m: " + num2str(tot_Orbcomm_frag_1m))
