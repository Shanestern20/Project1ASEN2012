clear;
clc;
close all;

% Load data (time,x,y)
dataA = readmatrix('Data_TCAS_A.csv'); 
dataB = readmatrix('Data_TCAS_B.csv');

% Separate data columns
timeA = dataA(:, 1);
xA = dataA(:, 2);
yA = dataA(:, 3);

timeB = dataB(:, 1);
xB = dataB(:, 2);
yB = dataB(:, 3);

%plot raw data positions
figure()
plot(xA,yA)
hold on
plot(xB,yB)




% xA and yA lin. fit
p_xA = polyfit(timeA, xA, 1); % p_xA(1) is the velocity uA, p_xA(2) is the xA0
p_yA = polyfit(timeA, yA, 1); % p_yA(1) is the velocity vA, p_yA(2) is the yA0

% xB and yB lin. fit
p_xB = polyfit(timeB, xB, 1); % p_xB(1) is the velocity uB, p_xB(2) is the xB0
p_yB = polyfit(timeB, yB, 1); % p_yB(1) is the velocity vB, p_yB(2) is the yB0



% get initial positions and velocities
xA0 = p_xA(2); yA0 = p_yA(2);
xB0 = p_xB(2); yB0 = p_yB(2);
uA = p_xA(1); vA = p_yA(1);
uB = p_xB(1); vB = p_yB(1);

%plot lines of best fit and extrapolate forward
atime = (min(timeA):700);
linea = zeros(length(atime),2);
for i = 1:length(atime)
    xApos = uA*atime(i) + xA0;
    yApos = vA*atime(i) + yA0;
    linea(i,1) = xApos;
    linea(i,2) = yApos;
end

plot(linea(:,1),linea(:,2))
btime = (min(timeB):700);
lineb = zeros(length(btime),2);
for i = 1:length(btime)
    xBpos = uB*btime(i) + xB0;
    yBpos = vB*btime(i) + yB0;
    lineb(i,1) = xBpos;
    lineb(i,2) = yBpos;
end
plot(lineb(:,1),lineb(:,2))

% Calculate tca
t_ca = (-( (xB0 - xA0) * (uB - uA) + (yB0 - yA0) * (vB - vA) ))/((uB - uA)^2 + (vB - vA)^2);


% Calculate min distance
xA_ca = xA0 + uA * t_ca;
yA_ca = yA0 + vA * t_ca;
xB_ca = xB0 + uB * t_ca;
yB_ca = yB0 + vB * t_ca;

%plot locations at tca
plot(xA_ca,yA_ca,'o')
plot(xB_ca,yB_ca,'o')

%min distance
D_min = sqrt((xB_ca - xA_ca)^2 + (yB_ca - yA_ca)^2);

%%%%%%%%%%%%%%%%finish error prop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assuming standard deviations (errors) for initial positions and velocities
sigma_xA0 = 0.05; % Example value in nmi
sigma_yA0 = 0.05; % Example value in nmi
sigma_xB0 = 0.05; % Example value in nmi
sigma_yB0 = 0.05; % Example value in nmi
sigma_uA = 0.01;  % Example value in nmi/s
sigma_vA = 0.01;  % Example value in nmi/s
sigma_uB = 0.01;  % Example value in nmi/s
sigma_vB = 0.01;  % Example value in nmi/s

% Partial derivatives of t_ca
partial_tca_xA0 = -(uB - uA) / ((uB - uA)^2 + (vB - vA)^2);
partial_tca_yA0 = -(vB - vA) / ((uB - uA)^2 + (vB - vA)^2);
partial_tca_xB0 = (uB - uA) / ((uB - uA)^2 + (vB - vA)^2);
partial_tca_yB0 = (vB - vA) / ((uB - uA)^2 + (vB - vA)^2);
partial_tca_uA = -((xB0 - xA0) * (2 * (uB - uA))) / ((uB - uA)^2 + (vB - vA)^2)^2;
partial_tca_vA = -((yB0 - yA0) * (2 * (vB - vA))) / ((uB - uA)^2 + (vB - vA)^2)^2;
partial_tca_uB = ((xB0 - xA0) * (2 * (uB - uA))) / ((uB - uA)^2 + (vB - vA)^2)^2;
partial_tca_vB = ((yB0 - yA0) * (2 * (vB - vA))) / ((uB - uA)^2 + (vB - vA)^2)^2;

% Error propagation for t_ca
sigma_tca = sqrt((partial_tca_xA0 * sigma_xA0)^2 + (partial_tca_yA0 * sigma_yA0)^2 + ...
                 (partial_tca_xB0 * sigma_xB0)^2 + (partial_tca_yB0 * sigma_yB0)^2 + ...
                 (partial_tca_uA * sigma_uA)^2 + (partial_tca_vA * sigma_vA)^2 + ...
                 (partial_tca_uB * sigma_uB)^2 + (partial_tca_vB * sigma_vB)^2);

% Error propagation for D_min (assume similar derivation steps)
partial_Dmin_xAca = (xA_ca - xB_ca) / D_min;
partial_Dmin_yAca = (yA_ca - yB_ca) / D_min;
partial_Dmin_xBca = -partial_Dmin_xAca;
partial_Dmin_yBca = -partial_Dmin_yAca;

sigma_Dmin = sqrt((partial_Dmin_xAca * sigma_xA0)^2 + (partial_Dmin_yAca * sigma_yA0)^2 + ...
                  (partial_Dmin_xBca * sigma_xB0)^2 + (partial_Dmin_yBca * sigma_yB0)^2);

% Display error results
fprintf('Error in t_ca: %.4f seconds\n', sigma_tca);
fprintf('Error in D_min: %.4f nmi\n', sigma_Dmin);

% Check for TCAS
TA_threshold = 3.3; 
RA_threshold = 2.0; 

if D_min < RA_threshold
    warning_msg = 'Resolution Advisory (RA) issued';
elseif D_min < TA_threshold
    warning_msg = 'Traffic Advisory (TA) issued';
else
    warning_msg = 'No advisory issued';
end

% Display Results
fprintf('Time of Closest Approach (t_ca): %.2f seconds\n', t_ca);
fprintf('Minimum Distance (D_min): %.2f nmi\n', D_min);
%fprintf('Error in t_ca......;
fprintf('Warning: %s\n', warning_msg);
