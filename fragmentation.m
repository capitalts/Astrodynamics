function [debris1]=fragmentation(ep, p1)

% reference: MASTER-8-Final-Report Section 2, namely 2.2.3.1 and 2.2.5
%               which revises the NASA SBM (EVOLVE 4.0) and other
%               previous studies

LB = 0.1;       % trackable diameter dtr = 10 cm; also, characteristic length L_c later

M_total_mass = p1.mass;

% FIX: Add all classes for objects
% cs=0.25*rand;
% reference: MASTER-8-Final-Report (pg 65/465, tab 2.1)
% Scaling factors for the revised NASA breakup model as used for future projections in EVOLVE 4.0
% Parent Object Type                                                      s
% PROTON (SL-12) ullage motors (SOZ units)                              0.1
% Other (non-SOZ) rocket bodies                                         1.0
% EORSATs (Soviet/Russian Electronic Ocean Reconnaissance Satellites)   0.6
% Molniya type early warning satellites                                 0.1
% All Soviet/Russian battery-related events                             0.5
% All Soviet/Russian anti-satellite tests (ASAT)                        0.3
% Other payloads                                                        1.0

% if (600<p1.mass) && (p1.mass<1000)
    cs = 1;      % scaling factor 's' in MASTER and NASA SBM
% else 
%     r = randi([1 6],1,1);      % randomly select type of original object
%     switch r
%         case 1
%             cs = 0.1;
%         case 2
%             cs = 1;
%         case 3
%             cs = 0.6;
%         case 4
%             cs = 0.1;
%         case 5 
%             cs = 0.5;
%         case 6
%             cs = 1;
%     end
% end

% reference: SPACE DEBRIS, Model and Risk Analysis
num = floor(6 * cs * LB ^(-1.6)); % eq(3.32);  LB is L_c

% reference: MASTER-8-Final-Report (pg 75/465, section 2.2.5.2)
% power law for explosions has been corrected and is now only used up to 1 m. 
% Above 1 m, only the randomly created large fragments are used.
global radiusearthkm
H = (p1.a*radiusearthkm)*(1-p1.ecco);
logH = log10(H);
if H <= 620 %[km]
    dtr = 8.9; %[cm]
elseif (620 < H) && (H <= 1300)
    dtr = 10^(-0.736748+0.604*logH); %[cm]
elseif (1300 < H) && (H <= 3800)
    dtr = 10^(-4.517+1.8186*logH); %[cm]
else % > 3800
    dtr = 1; %[m]
    dtr = dtr*100; %[cm]
end

% Henize factor (dtr in cm)
if dtr > 10^0.78 %[cm]
    fhz = sqrt( 10^ ( exp(-( (log10(dtr)-0.78)/0.637 )^2) ) ); % eq 2.8.2
else 
    fhz = sqrt(10);
end
num = num*fhz; % eq 2.8.1

% ORIGINAL METHOD for 'd' (transform to power law...?)
d=rand(floor(num+0.2*num),1);  %%%%% REVISIT;  sample diameter d uniformly..?
d=np_where(d < LB, 0*d, 1-d.^(-1.71)/ (LB ^(-1.71))); % remove d < LB, otherwise calculate diam; 
    %%%% ^^ REVISIT: exp refers to Collision dynamics..?
    % DEBUG ^
    %     figure(11); clf;
    %     plot(0:0.1:1, 1-[0:0.1:1].^(-1.71)/(LB^-1.71),'-x'); hold on;
    %     plot(0:0.1:1, 1-[0:0.1:1].^(-1.71)/((2*LB)^-1.71),'-x');

% NEW METHOD #1: uniformly sample in logspace (crude; produces multiple of same sized objects)
dd = logspace(log10(LB),log10(p1.radius*2),100);  % uniformly sample diameters in logspace
ndd = 6 * cs * dd.^-1.6;     % distribution of n's  (eq 2.68)
ndd = ndd / sum(ndd) * num;  % normalize to num from before
d = repelem(dd,ceil(ndd))';   % new d (for now)
%     figure(12); clf; dd=  0:0.1:10;
%     loglog(dd, 6*cs*dd.^-1.6,'-x'); hold on; % eq 2.68


% NEW METHOD #2: create PDF, then sample 'num' selections
dd = linspace(LB,p1.radius*2,1000);
ndd = 6 * cs * dd.^-1.6;     % distribution of n's  (eq 2.68)
d = repelem(dd,ceil(ndd))';   % new d (for now)


% NEW METHOD #3: do #2, but only up to 1m, then randomly sample larger objects
%                   as described in 2.2.3.1 "in case of an explosion event, the 
%                   power law diameter distribution is generally cut at 1 m and 
%                   a set of 4s to 8s remnant objects is randomly generated
%                   above this threshold." and in 2.2.5.2 

dd = linspace(LB,1,1000);   % sampled diameters ; power law valid from LB to 1 m
nddcdf =  dd.^-1.6;            % CUMULATIVE distribution of n's  (eq 2.68)
ndd = flip(diff(flip(nddcdf))); % diff to get PDF count
d_pdf = repelem(dd(1:end-1),round(ndd*10000))';   % PDF of debris objects between LB and 1 m
d = d_pdf(randperm(numel(d_pdf),round(num)));
    figure; subplot(411); plot(dd(1:end),nddcdf); xlabel('diam'); ylabel('CDF'); title('theoretical CDF (eq 2.68)');
    subplot(412); histogram(d_pdf); xlabel('diam'); ylabel('count'); title('PDF to sample from')
    subplot(413); semilogy(flip(cumsum(flip(histcounts(d_pdf,1000))))); xlabel('arb (diam)'); ylabel('cumulative count'); title('CDF of above');
    subplot(414); histogram(d); xlabel('diam'); ylabel('count'); title('PDF of sampled d');

    A = 0.556945*d.^(2.0047077);            % calculate area; Eq 2.72
    Am=zeros(length(d),1);
    dv=Am;    
    for k=1:length(d)
        Am(k,1) = func_Am(d(k), p1.objectclass);
        dv(k,1) = func_dv(Am(k))/1000;
    end    
    m = A./Am;                              % calculate mass
    m_remSum = M_total_mass - sum(m);       % remnant mass not accounted for by debris; divide this into 4 to 8 pieces (2.2.3.1)
    remDist = rand(randi([4,8]),1);
    m_rem = m_remSum * remDist / sum(remDist); % randomly distribute the remnant mass
    r_rem = (m_rem ./ p1.mass * p1.radius^3).^(1/3); % use same density (mass / r^3) as original object to get A_rem
    % test: (m_rem ./ r_rem.^3) / (p1.mass / p1.radius.^3)
    d_rem = r_rem * 2;
    d = [d; d_rem];
    A_rem = pi* r_rem.^2;           % need to calc these next few lines b/c d to m conversion invalid for the rems
    Am_rem = A_rem ./ m_rem;
    for k=1:length(Am_rem)
        dv_rem(k,1) = func_dv(Am_rem(k))/1000;
    end    


% EVOLVE4.0 method from C++ code
 % COMING SOON


d = d(d < 2*p1.radius); % filter by max size (diameter of original obj)
d = d(d > 1e-6);        % filter by min size

A = 0.556945*d.^(2.0047077);     % calculate area; Eq 2.72
Am=zeros(length(d),1);
dv=Am;

for k=1:length(d)
    Am(k,1) = func_Am(d(k), p1.objectclass);
    %dv(k,1) = func_dv(Am(k),'explosion')/1000;
    dv(k,1) = func_dv(Am(k))/1000;
end

m = A./Am;                      % calculate mass

% create samples on unit sphere
u=rand(length(dv),1)*2-1;
theta=rand(length(dv),1)*2*pi;

v = sqrt(1 - u.^2);
p = [v.*cos(theta) v.*sin(theta) u];

% dv velocity vectors
try
    dv_vec = [p(:,1).*dv  p(:,2).*dv  p(:,3).*dv];
catch
    return
end

fragments = [d A Am m dv dv_vec(:,1) dv_vec(:,2) dv_vec(:,3)];
% % DEBUG: HISTOGRAM OF THESE ^
%     figure(10);clf;
%     subplot(511); histogram(d); xlabel('d (m)');
%     subplot(512); histogram(A); xlabel('A (m^2)');
%     subplot(513); histogram(Am); xlabel('Am (m^2/kg)');
%     subplot(514); histogram(m); xlabel('m (kg)');
%     subplot(515); histogram(dv); xlabel('dv (m/s)');
%     subplot(511); title('fragment distributions', sprintf('Original mass: %0.1f kg, radius: %0.1f m',p1.mass, p1.radius));

idx = max(find(cumsum(m)<M_total_mass));    % choose index of last debris that cumsums to less than original mass

fragments = fragments(1:idx,:);

% distribute according to size first, then mass
var_foo = fragments(:,1) > p1.radius * 2;% fragments too large to originate from satelite 1
fragments1 = fragments(var_foo,:);
var_m1 = sum(fragments1(:,4)) ;
var_total_m = sum(fragments(:,4));

% if the second object is too small and there are no debris from it
if isempty(fragments(:,1) < p1.radius * 2)
%     fragments2 = [];
else
    frag1=fragments((fragments(:,1) < p1.radius * 2),:);
    r = randperm(size(frag1,1));
    var_rest = frag1(r,:);
    var_cumsum = cumsum(var_rest(:,4)); % calculate the cumulative sum of the mass of the rest
    var_idx = max(find(var_cumsum > (p1.mass/(M_total_mass) * var_total_m) - var_m1));
    fragments1 = [fragments1 ; var_rest(1:var_idx,:)];
%     fragments2 = var_rest(var_idx:end,:);
end

% DEBUG: HISTOGRAM OF THESE ^
%     figure(10);clf;
%     subplot(511); fragments1(:,1); xlabel('d (m)');
%     subplot(512); fragments1(:,2); xlabel('A (m^2)');
%     subplot(513); fragments1(:,3); xlabel('Am (m^2/kg)');
%     subplot(514); fragments1(:,4); xlabel('m (kg)');
%     subplot(515); fragments1(:,5); xlabel('dv (m/s)');
%     subplot(511); title('fragment distributions', sprintf('Original mass: %0.1f kg, radius: %0.1f m',p1.mass, p1.radius));


%debris1 = func_create_tles(ep, p1, fragments1);
debris1 = func_create_tlesv2(ep, p1, fragments1);

[debris1]=add_mass_radius(debris1);

% Initialize
for k=1:length(debris1)
    debris1{k}.controlled = 0;
    debris1{k}.constel = 0;
    debris1{k}.objectclass='Rocket Fragmentation Debris';%FIX: Should depend on which object fragments
end


    
