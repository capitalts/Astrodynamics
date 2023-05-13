%% Jackie H Smith and Tory D Smith, MIT
% Astrodynamics, Class Project
% 16 May 2023

close all; clear; clc

%% Russian ASAT Test
% Date: Nov 15, 2021
% Mass (COSMOS 1408): 1750 kg
% Velocity (ASAT): 27,000 kilometers per hour
% Relative Velocity: 4.6 km/s (3-6 km/s)
% Source: https://www.space.com/russia-anti-satellite-missile-test-first-of-its-kind
% https://www.armscontrol.org/act/2021-12/news/russian-asat-test-creates-massive-debris
% https://www.agi.com/blog/2021/12/asat-weapon-interception-debris-field

mass_secondary = [50, 63.7, 100]; % kg
mass_primary = 1750; % kg
velocity_relative = [1*1000, 4.6*1000, 6*1000]; % m/s

% energy equation
% https://github.com/nasa/CARA_Analysis_Tools/blob/master/conjunction_consequence_assessment/References/Hejduk%202017.pdf

figure(1)
for j = 1:length(velocity_relative)
    for i = 1:length(mass_secondary)
        E = mass_secondary(i)*velocity_relative(j)^2/(2*mass_primary); 
        
        if E < 40000
            % non-catastrophic
            P = mass_secondary(i)*velocity_relative(j)/1000; % km/s
            disp(['The collision is non-catastrophic at ASAT mass ' num2str(mass_secondary(i)) ' kg and rel_vel ' num2str(velocity_relative(j)) ' m/s.']);
        else
            % catastrophic
            P = mass_secondary(i)+mass_primary;
            disp(['The collision is catastrophic at ASAT mass ' num2str(mass_secondary(i)) ' kg and rel_vel ' num2str(velocity_relative(j)) ' m/s.']);
        end
        
        Lc = 0.01:0.01:10;
        
        N = colCons(Lc', P);
        
        loglog(Lc, N, 'Linewidth', 2)
        hold on
    
    end
end
title('Debris pieces above Characteristic Length', 'Fontsize', 24)
xlabel('Characteristic Length (m)', 'Fontsize', 16)
ylabel('# of debris larger than L_c', 'Fontsize', 16)
legend('ASAT 50 kg - vel 1km/s', 'ASAT 63.7kg - vel 1km/s', 'ASAT 100kg - vel 1km/s',...
    'ASAT 50 kg - vel 4.6km/s', 'ASAT 63.7kg - vel 4.6km/s', 'ASAT 100kg - vel 4.6km/s',...
    'ASAT 50 kg - vel 6km/s', 'ASAT 63.7kg - vel 6km/s', 'ASAT 100kg - vel 6km/s')
ylim([1 10^5])
grid on


%% Functions
function N = colCons(Lc, P)
    N = 0.1*P^0.75*Lc.^-1.71;
end