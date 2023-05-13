function [debrisMR]=frag_exp_SBM_C(ep, p1)

% Modified from fragmentation.m by djang in Oct 2022

% Fragmentation (Explosion) model following 
% NASA EVOLVE 4.0 standard breakup model (2001)
% with the revision in ODQN "Proper Implementation of the 1998 NASA Breakup
% Model" (2011)

% Explosions in the SBM follows the power law of the independent variable:
% the characteristic length (Lc).  Mass is controlled by adding in large
% pieces of debris that account for most of the mass of the original object

% Quote: 
% In the region above 1 m the data shows the deposit of several large pieces
% that do not necessarily follow the power law distribution. These would realistically be larger,
% more massive components farther from the explosion center (e.g., remnants of equipment
% shelves, pressurant tanks, nozzle bells, etc.). In fact these fragments account for the bulk of
% fragment mass. Based on this understanding, the correct implementation of the 1998 NASA
% Breakup Model includes the distribution of fragments from 1 mm to 1 m following the
% power law distribution in Equation 1, with an additional two to eight large fragments after
% 1 m, keeping mass conserved

% Note that the SBM distribution is modified by ESA for MASTER-2009 and MASTER-8
% and also by other authors in other papers e.g. Cimmino (2021), etc

% Plotting the power-law distribution of explosions by SBM; compare to Fig 1 of ODQN
% In the SBM, this region is the same for all explosions between 0.1 and 1 m
%     Lcs = logspace(-2,1,1000);
%     loglog(Lcs, 6 * Lcs.^-1.6); grid on
%     ylim([1,10000]);  set(gca, 'XTickLabel',get(gca,'XTick')) ; set(gca, 'YTickLabel',get(gca,'YTick')) 
%     legend('NASA: Ncum = 6 Lc^-1.6');

LB = 0.1;       % trackable diameter dtr = 10 cm; also, characteristic length L_c later

% M_total_mass = p1.mass;
% 
% % FIX: Add all classes for objects
% % cs=0.25*rand;
% % reference: MASTER-8-Final-Report (pg 65/465, tab 2.1)
% % Scaling factors for the revised NASA breakup model as used for future projections in EVOLVE 4.0
% % Parent Object Type                                                      s
% % PROTON (SL-12) ullage motors (SOZ units)                              0.1
% % Other (non-SOZ) rocket bodies                                         1.0
% % EORSATs (Soviet/Russian Electronic Ocean Reconnaissance Satellites)   0.6
% % Molniya type early warning satellites                                 0.1
% % All Soviet/Russian battery-related events                             0.5
% % All Soviet/Russian anti-satellite tests (ASAT)                        0.3
% % Other payloads                                                        1.0
% 
% % if (600<p1.mass) && (p1.mass<1000)
%     cs = 1;      % scaling factor 's' in MASTER and NASA SBM
% % else 
% %     r = randi([1 6],1,1);      % randomly select type of original object
% %     switch r
% %         case 1
% %             cs = 0.1;
% %         case 2
% %             cs = 1;
% %         case 3
% %             cs = 0.6;
% %         case 4
% %             cs = 0.1;
% %         case 5 
% %             cs = 0.5;
% %         case 6
% %             cs = 1;
% %     end
% % end
% 
% % reference: SPACE DEBRIS, Model and Risk Analysis
% % Cumulative number of debris between LB and 1 m diameter
% num = floor(6 * cs * LB ^(-1.6));       % eq(3.32);  LB is L_c;  
% 
% 
% % NEW METHOD #1: uniformly sample in logspace (crude; produces multiple of same sized objects)
% % dd = logspace(log10(LB),log10(p1.radius*2),100);  % uniformly sample diameters in logspace
% % ndd = 6 * cs * dd.^-1.6;     % distribution of n's  (eq 2.68)
% % ndd = ndd / sum(ndd) * num;  % normalize to num from before
% % d = repelem(dd,ceil(ndd))';   % new d (for now)
% % %     figure(12); clf; dd=  0:0.1:10;
% % %     loglog(dd, 6*cs*dd.^-1.6,'-x'); hold on; % eq 2.68
% 
% % NEW METHOD #2: create PDF, then sample 'num' selections
% % dd = linspace(LB,p1.radius*2,1000);
% % ndd = 6 * cs * dd.^-1.6;     % distribution of n's  (eq 2.68)
% % d = repelem(dd,ceil(ndd))';   % new d (for now)
% 
% 
% % NEW METHOD #3: do #2, but only up to 1m, then randomly sample larger objects
% %                   as described in 2.2.3.1 "in case of an explosion event, the 
% %                   power law diameter distribution is generally cut at 1 m and 
% %                   a set of 4s to 8s remnant objects is randomly generated
% %                   above this threshold." and in 2.2.5.2 
% 
% dd = linspace(LB,1,1000);   % sampled diameters ; power law valid from LB to 1 m
% nddcdf = 6 * cs * dd.^-1.6;            % CUMULATIVE distribution of n's  (eq 2.68)
% ndd = flip(diff(flip(nddcdf))); % diff to get PDF count
% d_pdf = repelem(dd(1:end-1),ceil(ndd*10000))';   % PDF of debris objects between LB and 1 m
% d = d_pdf(randperm(numel(d_pdf),ceil(num)));
% % note that this dist is the same for all objects no matter mass or size of parent object
% 
% %     figure(22); clf; subplot(321); 
% %     Lcs = logspace(-2,1,100); loglog(Lcs, 6 * Lcs.^-1.6); grid on; ylim([1,10000]); 
% %     hold on; loglog(dd(1:end),nddcdf); xlabel('diam'); ylabel('CDF'); title('theoretical reverse CDF');
% %     legend('NASA: N_{cum} = 6s Lc^{-1.6}','nddcdf from code');
% %     subplot(323); histogram(d_pdf); xlabel('diam'); ylabel('count'); title('PDF to sample from')
% %     subplot(325); loglog(Lcs(1:end-1), flip(cumsum(flip(histcounts(d_pdf,Lcs)))),'-x'); xlabel('Diam (m)'); ylabel('cumulative count'); title('CDF of above (for shape)');
% %     subplot(324); histogram(d); xlabel('diam'); ylabel('count'); title('PDF of sampled d');
% %     subplot(326); loglog(Lcs, 6 * Lcs.^-1.6); hold on; loglog(Lcs(1:end-1), flip(cumsum(flip(histcounts(d,Lcs)))),'-x'); xlabel('Diam (m)'); ylabel('cumulative count'); title('CDF of sampled debris');
%     
% % calculate mass of objects [LB, 1 m] by d > A > Am > m
% A = 0.556945*d.^(2.0047077);            % calculate area; Eq 2.72
% Am = func_Am(d, p1.objectclass);
% m = A./Am;
% % sum(m);       % about 100 kg is accounted for debris < 1 m
% 
% 
% if sum(m) < M_total_mass        % if debris mass is less than total mass,
%     % Randomly select "additional two to eight large fragments after 1 m,
%     % keeping mass conserved" per ODQN.  "these fragments account for the bulk
%     % of fragment mass"
% 
%     % Assign remnant mass > d > A > Am > dv
%     m_remSum = M_total_mass - sum(m);       % remnant mass
%     remDist = rand(randi([2,8]),1);         % randomly distribute by random number [2,8] <<<<<<<
%     m_rem = m_remSum * remDist / sum(remDist); % randomly distribute the remnant mass
%     d_rem = (m_rem ./ p1.mass * p1.radius^3).^(1/3) * 2; % use same density (mass / r^3) as original object to get d_rem <<<<<<<<<<
%     % test: (m_rem ./ r_rem.^3) / (p1.mass / p1.radius.^3)
%     A_rem = 0.556945*d_rem.^(2.0047077);     % calculate area; Eq 2.72
%     Am_rem = A_rem ./ m_rem;
% else
% %     % method 1) sort by mass, keep those until mass adds up (bottom-up) <<<<<<<<<<<<<<<
% %     [msort, mord] = sort(m);
% %     lastidx = find(cumsum(msort) < M_total_mass,1,'last');  % last valid index of sorted list
% %     valididx = mord(1:lastidx);
% %     m = m(valididx);
% %     d = d(valididx);
% %     A = A(valididx);
% %     Am = Am(valididx);
%     % method 2) sort by Lc, keep smallest objects until mass adds up <<<<<<<<<<<<<<
%     [dsort, dord] = sort(d);
%     lastidx = find(cumsum(m(dord)) < M_total_mass,1,'last');  % last valid index of sorted list
%     valididx = dord(1:lastidx);
%     m = m(valididx);
%     d = d(valididx);
%     A = A(valididx);
%     Am = Am(valididx);
%     
%     d_rem = []; A_rem = []; Am_rem = []; m_rem = [];  % no "remnants" exist
% end
% 
% 
% % NASA SBM method directly from C++ code: COMING SOON <<<<<<<<<<<<<
% 
% 
% % Assign dv to random directions; create samples on unit sphere
% dv = func_dv([Am; Am_rem])/1000;       % km/s
% u=rand(length(dv),1)*2-1;
% theta=rand(length(dv),1)*2*pi;
% 
% v = sqrt(1 - u.^2);
% p = [v.*cos(theta) v.*sin(theta) u];
% 
% % dv velocity vectors
% try
%     dv_vec = [p(:,1).*dv  p(:,2).*dv  p(:,3).*dv];
% catch
%     warning('error assigning dv to unit sphere')
%     return
% end


% % create fragments
% fragments = [[d; d_rem] [A; A_rem] [Am; Am_rem] [m; m_rem] dv dv_vec(:,1) dv_vec(:,2) dv_vec(:,3)];
% % DEBUG: HISTOGRAM OF THESE ^
% %     figure(10);clf;
% %     subplot(511); histogram(d,100); xlabel('d (m)');
% %     subplot(512); histogram(A,100); xlabel('A (m^2)');
% %     subplot(513); histogram(Am,100); xlabel('Am (m^2/kg)');
% %     subplot(514); histogram(m,100); xlabel('m (kg)');
% %     subplot(515); histogram(dv * 1000,[0:10:1000]); xlabel('dv (m/s)');
% %     subplot(511); title('fragment distributions', sprintf('Original mass: %0.1f kg, radius: %0.1f m',p1.mass, p1.radius));
% %     pause;
% 
% % below: no need, as total mass should equal M_total_mass
% % idx = max(find(cumsum(m)<M_total_mass));    % choose index of last debris that cumsums to less than original mass
% % fragments = fragments(1:idx,:);
% if sum(m)>M_total_mass
%     warning('Sum of fragment mass (%0.1f kg) exceeds total mass of original object (%0.1f kg)', ...
%         sum(m), M_total_mass);
% end
% 
% if any(fragments(:,1) > p1.radius*2)
%     warning('Some fragments exceed the diameter of the original object (%0.1f m)', ...
%         p1.radius*2);
% end

% % % distribute according to size first, then mass
% % var_foo = fragments(:,1) > p1.radius * 2;% fragments too large to originate from satelite 1
% % fragments1 = fragments(var_foo,:);
% % var_m1 = sum(fragments1(:,4)) ;
% % var_total_m = sum(fragments(:,4));
% % 
% % % if the second object is too small and there are no debris from it
% % if isempty(fragments(:,1) < p1.radius * 2)
% % %     fragments2 = [];
% % else
% %     frag1=fragments((fragments(:,1) < p1.radius * 2),:);
% %     r = randperm(size(frag1,1));
% %     var_rest = frag1(r,:);
% %     var_cumsum = cumsum(var_rest(:,4)); % calculate the cumulative sum of the mass of the rest
% %     var_idx = max(find(var_cumsum > (p1.mass/(M_total_mass) * var_total_m) - var_m1));
% %     fragments1 = [fragments1 ; var_rest(1:var_idx,:)];
% % %     fragments2 = var_rest(var_idx:end,:);
% % end

[fragments] = Cversion_NASAbreakup(LB,'EXPLOSION',p1.mass,p1.v,p1.objectclass);

%debris1 = func_create_tles(ep, p1, fragments1);
debris1 = func_create_tlesv2(ep, p1, fragments);   % create orbits (OEs) and add mass and radius

% [debrisMR] = add_mass_radius(debris1);
%%% ^^^^^^^^^ need to compare this to mass and radius produced from within this script
% note: original add_mass_radius doesn't look right
% todo: DEBUG TEST: COMPARE: fragments, debris1 (input) , debrisMR (output)

debrisMR = debris1;


% Initialize fields
frag_objectclass1 = filter_objclass_fragments(p1.objectclass); %Assign object class to fragments according to parent particle
for k=1:length(debrisMR)
    debrisMR{k}.controlled = 0;
    debrisMR{k}.constel = 0;
    debrisMR{k}.objectclass=frag_objectclass1;
end


    
