clear all
close all
clc
config_name = 'optimization';

DIRECTORY = '..\optimization\';

%% Mirror array params
% N = 401; % mirror rows - for now, enfore Odd
% M = 401;% mirror cols - for now, enfore Odd
N = 731; % mirror rows - for now, enfore Odd
M = 731;% mirror cols - for now, enfore Odd
% % rotation params - not yet implemented
% phi_ipd = 0; % in plane rotation (degrees)


%% Last modified - PR - Feb 7,2021

%% Hard limits for X & Y axes when displaying PSF & PSF Slice
XAXIS_LIMITS = 40; % 2mm
YAXIS_LIMITS = 40; % 2mm

%% Load MAT file
% for frac_order_p = 0:0.1:0.5
%     for frac_order_q = -0.5:0.1:0.5
frac_order_p = 0.3;%0.2;  % Specified in fractional order (-1,1) for MWIR
frac_order_q = 0.3;
strBaseFileName = sprintf('BeamSteer_p=%1.2f,q=%1.2f',frac_order_p,frac_order_q);
% Replace decimal point with 'd'
strBaseFileName = replace(strBaseFileName,'.','d');
%
load([DIRECTORY strBaseFileName '.mat']);

%% physical mirror params
PITCH = 13.68e-6; % mirror pitch (m)
FF = 0.92; % mirror fill factor
% R_m = @(lambda) 0.96; % mirror refelctivitiy function
R_m = 0.96*ones(N,M);
% mirror height (mirror diagonal length/sqrt(2) * sind(12))
H = data.mirror_height * ones(N,M); % loaded from MAT file
theta_dist = 'gaussian'; % 
theta_nom = 12; % nominal tilt angle 
theta_std = 0; % tilt angle stddev

%% substrate params
R_s = @(x) 0.6; % Substrate refelctivitiy function

%% Dervived params
A_s = N*M*PITCH;
% build mirror dist
switch theta_dist
  case 'gaussian'
    m_theta = (theta_std^2*(randn(N,M))+theta_nom);
  otherwise
    error('Unrecognized tilt angle distrubition')
end

%% Illumination Params
lambda = 4e-6; % wavelength (m)
u = data.beta_i;                                  % u = beta
v = sqrt(1-(data.alpha_i.^2+data.beta_i.^2));     % v = gamma
azel = uv2azel([u;v]);
az_i = azel(1);  % x to y axis 
el_i = azel(2);  % y to z axis
fprintf('Incident Illumination - Azimuth,Elevation: %2.3f,%2.3f\n',az_i,el_i);
u = data.beta_tgt;                                % u = beta
v = sqrt(1-(data.alpha_tgt.^2+data.beta_tgt.^2)); % v = gamma
azel = uv2azel([u;v]);
az_tgt = azel(1);  % x to y axis 
el_tgt = azel(2);  % y to z axis
fprintf('Target - Azimuth,Elevation: %2.3f,%2.3f\n',az_tgt,el_tgt);

bitmask = data.bitmask(1:M,1:N);
% bitmask = 0*bitmask;
% bitmask = data.bitmask(round(M/2),round(N/2));
figure('Name','Optimized Result - BITMASK'), 
  imshow(bitmask,'InitialMagnification','fit'); 
  colorbar; colormap gray;
  caxis([-1 1]);
  set(gca,'YDir','normal');
  title('Bit Mask DMD');
  
%% build param structs
% mirror array params
pparams = struct(...
    'pitch',PITCH,...
    'ff',FF,...
    'R_m',R_m,...
    'R_s',R_s,...
    'N',N,...
    'M',M,...
    'h',0*H,...
    'm_theta',m_theta,...
    'ip_rot',0 ...
);
n_start = ((-(N-1)/2));
m_start = ((-(M-1)/2));
[mirror_n, mirror_m] = meshgrid((1:N)+n_start-1,(1:M)+m_start-1);
center_n = 0;
center_m = 0;
beam_waist_n = 365.4971/2;%N/2; % in mirrors
beam_waist_m = 365.4971/2;%M/2; % in mirrors

R = [cosd(az_i) -sind(az_i); sind(az_i) cosd(az_i)];
mirrors = R*[(mirror_n(:)-center_n)';(mirror_m(:)-center_m)'];
mirrors_n = reshape(mirrors(1,:),N,M);
mirrors_m = reshape(mirrors(2,:),N,M);

% E_m = exp(-pi*(((mirror_n-center_n)/(beam_waist_n)).^2 + ((mirror_m-center_m)/(beam_waist_m)).^2));
E_m = ...
  exp( ...
    -pi*( ...
      (mirrors_n/beam_waist_n).^2 + ...
      ((mirrors_m/beam_waist_m)/sind(el_i)).^2));
%    
figure('Name','Illumination Intensity'), 
  imshow(E_m,'InitialMagnification','fit' );
set(gcf,'Units','Normalized');
figPos = get(gcf,'OuterPosition');
  
%
E_src = sum(E_m(:))* A_s;
I_src = E_src.^2 ;

% %%
% az_in = 45;       % x to y axis 
% el_in = (90-24);  % y to z axis
% uv_in = azel2uv([az_in;el_in]);
% % % dervived params
% beta_i = uv_in(1); 
% gamma_i = uv_in(2);
% alpha_i = sqrt(1-(beta_i^2+gamma_i^2));
% %%

% illumination params
illumination_params = struct(...
    'lambda',lambda,...
    'alpha_i',data.alpha_i,...data.
    'beta_i',data.beta_i,...
    'center_n',center_n,...
    'center_m',center_m,...
    'beam_waist_n', beam_waist_n,...
    'beam_waist_m', beam_waist_m,... 
    'E_m',E_m...
);

%% View parameters
r = data.R_0;       % Distance @ which to create spot
% r = 2.5/100;       % Distance @ which to create spot
%{
view_params = struct(...
    'view_style','full',... % 'abs','order','full'
    'r',r,...
    'beta_min',data.beta_tgt - 2.5e-3,...
    'beta_max',data.beta_tgt + 2.5e-3,...
    'alpha_min',data.alpha_tgt - 2.5e-3,... 
    'alpha_max',data.alpha_tgt + 2.5e-3,...
    'order_x',0,...
    'order_y',0,...
    'spacing_param',0.05,...
    'units','dcs',... % 'dcs','rad','deg'
    'samples',512 ... prefer power of 2, max 1024 right now, need to subdivide grid in code for handling denser sampline since it is tied to cuda threadblock size
);
%}
view_params = struct(...
    'view_style','order',... % 'abs','order','full'
    'r',r,...
    'order_x',data.frac_order_p,...
    'order_y',data.frac_order_q,...
    'spacing_param',0.025,...125
    'units','dcs',... % 'dcs','rad','deg'
    'samples',512 ... prefer power of 2, max 1024 right now, need to subdivide grid in code for handling denser sampline since it is tied to cuda threadblock size
);

%% calculate field and intensity pattern
tic
[dalpha,dbeta,U] = ...
  dmd_diffraction(pparams,illumination_params,bitmask,view_params,false);
toc

% [dalpha,dbeta,L] = ...
%   dmd_diffracted_2(pparams,illumination_params,bitmask,view_params,false);
% toc

%{
% clear all
load('phasors.mat')
input('enter to begin plotting')
for i = 1:size(phasors,3)
    s = phasors(:,:,i);
    figure(10), scatter(real(s(:)),imag(s(:))), title(['mirror: ',num2str(i)]), xlabel('Re(U)'), ylabel('Im(U)');
    pause(0.1)
end
%}

%% Diffracted intensity
I = abs(U).^2;

%% irridance
[Alpha,Beta] = meshgrid(dalpha,dbeta);
Gamma = sqrt(1-Alpha.^2 - Beta.^2);
g_real = imag(Gamma)==0;
% L_tot = trapz(dbeta,trapz(dalpha,double(L),2));
% L_real = trapz(dbeta,trapz(dalpha,double(L).*g_real,2));
% K = L_tot/L_real;
% U_tot = double(sum(abs(L(:)).^2));
% U_real = double(sum(abs(L(g_real)).^2));
% K = U_tot/U_real
% L_p = K*L.*g_real; %*gamma_i*(illumination_params.lambda^2)/A_s
% L_p_tot = trapz(dbeta,trapz(dalpha,double(L_p),2));

% max_L = max(L_p(:));

% I_L = L_p.*Gamma.*A_s;


%% Last modified - PR - Feb 7,2021

%% Compute spot-size, eccentricity & centroid
% insensitive to diffracted energy normlization
[nR,nC] = size(I); 
% Inter-sample-spacing in physical coordinates
ISS = 2 * view_params.spacing_param ./ view_params.samples .* view_params.r; 
% Fit ellipse to diffracted intensity
In = I/max(I(:));
% In = double(nthroot(In,4));
% th = graythresh(In);
th = 0.01;  % 20 dB below peak-value
bw = imbinarize(In,th);
stats = regionprops(bw,'all');
[~,idx]=max(arrayfun(@(x) x.Area, stats));
stats = stats(idx);
% Eccentricity
% fprintf('Eccentricity = %f\n',stats.Eccentricity);
% Centroid
xbar = ISS * (stats.Centroid(1)-((nC-1)/2));
ybar = ISS * (stats.Centroid(2)-((nR-1)/2));
% Ellipse params
a = ISS * stats.MajorAxisLength/2;
b = ISS * stats.MinorAxisLength/2;
%
SpotSize_Major = ISS * stats.MajorAxisLength;
SpotSize_Minor = ISS * stats.MinorAxisLength;
%
theta = pi*stats.Orientation/180;
R = [ cos(theta)   sin(theta)
      -sin(theta)   cos(theta)];  % flipped because of image coordinates
%
phi = linspace(0,2*pi,50);
%
xy = [a*cos(phi); b*sin(phi)];
xy = R*xy;
xe = xy(1,:) + xbar;
ye = xy(2,:) + ybar;
%
xg = -(nC-1)/2:(nC-1)/2;  xg = xg * ISS;
yg = -(nR-1)/2:(nR-1)/2;  yg = yg * ISS;
%
mm = 1e-3;
figure(102);
set(gcf,'units','Normalized','outerposition',[0 0 1 1]);
%   imshow(nthroot(In,2),[],'InitialMagnification','fit','XData',xg/mm,'YData',yg/mm);
  imagesc(nthroot(In,1),'XData',xg/mm,'YData',yg/mm); colorbar, colormap gray
  set(gca,'YDir','normal');
%   hold on;
%   scatter(data.alpha_tgt - data.alpha_i, data.beta_tgt-data.beta_i,'r')
  hold on
  
  line([xbar,xbar]/mm,[yg(1),yg(end)]/mm,'Color',[255,165,0]/255,'LineWidth',2);
  line([xg(1),xg(end)]/mm,[ybar,ybar]/mm,'Color',[255,165,0]/255,'LineWidth',2);
  plot(xe/mm,ye/mm,'r','LineWidth',2);
%   plot(xbar/mm,ybar/mm,'.','markersize',16);
  prev_xax_limits = xlim();  prev_yax_limits = ylim();
  % Reset x & y axis limits to enable cropping
  xlim([-1,1]*XAXIS_LIMITS); ylim([-1,1]*YAXIS_LIMITS);
  %
  title('Diffracted Radiance','fontsize',18); xlabel('millimeters'); ylabel('millimeters');
  set(gca,'fontname','Candara','fontsize',14);
  % Print figure window
  print('-dpng','-r300',[DIRECTORY strBaseFileName '_PSF.png']);
  % Trim the image by discrading outrer rows & columns of pixels with constant intensity
  boolOverwrite = true;
  [TC,TR] = TrimImage([DIRECTORY strBaseFileName '_PSF.png'], boolOverwrite);
  %{
  subplot(122);
    imshow(angle(U)/pi,[],'InitialMagnification','fit','XData',xg,'YData',yg);  
    title('Phase')
  %}
  % Reset x & y axis limits to enable cropping
  xlim(prev_xax_limits); ylim(prev_yax_limits);
set(gcf,'OuterPosition',figPos);  % Reset figure to standard size

% error

%% Compute peak to side lobe ratio from diffracted intensity
% insensitive to diffracted energy normlization
I = abs(U).^2;    
In = I/max(I(:));
[nR,nC] = size(I);
I_beta_slice = mean(In(:,[view_params.samples/2,view_params.samples/2+1]),2);
I_alpha_slice = mean(In([view_params.samples/2,view_params.samples/2+1],:),1);
%
% Compute peak to side-lobe ratio in X direction
It = 10*log10(I_alpha_slice);
[mx,mxloc] = max(It);
[pks,locs] = findpeaks(It);
[~,ord] = sort(abs(pks-mx),'ascend');
sidelobe_loc = ord(2);
PSL_X = abs(pks(sidelobe_loc)-mx); % expressed in dB
%
% Compute peak to side-lobe ratio in X direction
It = 10*log10(I_beta_slice);
[mx,mxloc] = max(It);
[pks,locs] = findpeaks(It);
[~,ord] = sort(abs(pks-mx),'ascend');
sidelobe_loc = ord(2);
PSL_Y = abs(pks(sidelobe_loc)-mx); % expressed in dB
%
figure,
set(gcf,'units','Normalized','outerposition',[0 0 1 1]);
  semilogy(xg/mm, I_alpha_slice,'-','LineWidth',4); 
  prev_xax_limits = xlim();
  xlim([-1,1]*XAXIS_LIMITS);   % Reset x limit to enable cropping 
%   xlabel('millimeters');
    set(gca,'XTickLabel',[])
  set(gca,'fontname','Candara','fontsize',96,'color','none');
%   daspect([1,0.75,1]);
  % Print figure window
  print('-dpng','-r300',[DIRECTORY strBaseFileName '_XSlice.png']);
  % Trim the image by discrading outrer rows & columns of pixels with constant intensity
  boolOverwrite = true;
  [TC,TR] = TrimImage([DIRECTORY strBaseFileName '_XSlice.png'], boolOverwrite);
  xlim(prev_xax_limits);
  daspect([1,0.25,1]);
set(gcf,'OuterPosition',figPos);  % Reset figure to standard size
%
figure, 
set(gcf,'units','Normalized','outerposition',[0 0 1 1]);
  semilogy(yg/mm, I_beta_slice,'-','LineWidth',4); 
  prev_xax_limits = xlim();
  xlim([-1,1]*XAXIS_LIMITS);   % Reset x limit to enable cropping 
  xlabel('millimeters');
  set(gca,'fontname','Candara','fontsize',36,'color','none');
%   daspect([1,0.75,1]);
  % Print figure window
  print('-dpng','-r300',[DIRECTORY strBaseFileName '_YSlice.png']);
  % Trim the image by discrading outrer rows & columns of pixels with constant intensity
  boolOverwrite = true;
  [TC,TR] = TrimImage([DIRECTORY strBaseFileName '_YSlice.png'], boolOverwrite);
  xlim(prev_xax_limits);
  daspect([1,0.25,1]);
set(gcf,'OuterPosition',figPos);  % Reset figure to standard size
  
%% Compute throughput
% fprintf('Throughput %f\n',NaN);

%% Write DMD pattern to disk
imwrite(0.5*(1+bitmask),[DIRECTORY strBaseFileName '_DMD_Pattern.png']); 

%% Write to excel file
Az_Incident = az_i;
El_Incident = el_i;
Az_Target = az_tgt;
El_Target = el_tgt;
%
Eccentricity = stats.Eccentricity;
Throughput = NaN;                       % Danyal TODO
Peak_to_SideLobe_Ratio = min(PSL_X,PSL_Y);
%
DMD_BitMask = {[DIRECTORY strBaseFileName '_DMD_Pattern.png']};
PSF = {[DIRECTORY strBaseFileName '_PSF.png']};
MAT_File = {[DIRECTORY strBaseFileName '.mat']};


sim_res = struct('physical_params',pparams,...
              'illumination_params',illumination_params,...
              'bitmask',bitmask,...
              'view_params',view_params,...
              'alpha_axis',dalpha,...
              'beta_axis',dbeta,...
              'I_U',I,...
              'U',U);
%               'L',L,...
%               'L_p',L_p,...

%               'I_L',I_L,...f
              

save(MAT_File{1},'sim_res','-append');

%
Tbl = ...
  table(...
    Az_Incident,El_Incident,...
    Az_Target,El_Target,...
    Throughput,Eccentricity,...
    Peak_to_SideLobe_Ratio,...
    SpotSize_Major,...
    SpotSize_Minor,...
    DMD_BitMask,...
    PSF,...
    MAT_File...
    );
% Write table to Excel
filename = 'NAVAIR_SBIR_Perf_Metrics.xlsx';
writetable(Tbl,[DIRECTORY filename],'WriteMode','append','Sheet',1);

close all
% clear all

%     end
% end

%% Save data back to MAT file - TODO
% Use save(,'-append');

error('DONE')

%% Unused









figure(101), 
%   imshow(nthroot(I,4),[],'XData',az(1,:),'YData',el(:,1),'InitialMagnification','fit');
%   imshow(nthroot(I,4),[],'InitialMagnification','fit');
  imshow(nthroot(I,4),[],'XData',dalpha,'YData',dbeta,'InitialMagnification','fit');
  xlabel('\alpha','fontsize',18)
  ylabel('\beta','fontsize',18)
  set(gca,'YDir','normal');
  axis square
  colorbar

% Fitting a surface to Diffracted radiance
In = double(nthroot(I,4)); In = In/max(In(:));
[Alpha,Beta]= meshgrid(dalpha,dbeta);

sf = fit([Alpha(:), Beta(:)], In(:), 'poly22'); 


error('DONE')
  
  
  
%% Visualization
%{
% Draw hemisphere wireframe
[Alpha,Beta]= meshgrid(dalpha,dbeta);
% Modified by PR - capitalized Gamma
Gamma = sqrt(1-(Alpha.^2 + Beta.^2));
greal = (imag(Gamma)==0);
Gamma = Gamma .* greal;
%
az_el = uv2azel([Beta(:)';Gamma(:)']);
az = az_el(1,:);
el = az_el(2,:);
maxAzimuth = max(az);     minAzimuth = min(az);
maxElevation = max(el);   minElevation = min(el);
%
[X,Y,Z] = sphere(500);
el = atan2d(Z,sqrt(X.^2+Y.^2));
az = atan2d(Y, X);
ind = (az >= minAzimuth & az <= maxAzimuth) & ...
      (el >= minElevation & el <= maxElevation); % // Find those indices
X2 = X; Y2 = Y; Z2 = Z; %// Make a copy of the sphere co-ordinates
X2(~ind) = NaN; 
Y2(~ind) = NaN; 
Z2(~ind) = NaN; %// Set those out of bounds to NaN
%
surf(r*X2,r*Y2,r*Z2,'FaceColor','white','FaceAlpha',0.5); %// Base sphere;
hold on;
%}

% [X,Y,Z] = sphere(20);
% r = data.R_0;
% surf(r*X,r*Y,r*Z,'FaceColor','white','FaceAlpha',0.5); %// Base sphere;
% hold on;

% Overlay spot on sphere
In = nthroot(I,4); In = In/max(In(:));
r = data.R_0 + 0.0*In;
X = r .* Alpha;
Y = r .* Beta;
Z = r .* Gamma;
% surf(X,Y,Z);
surf(X,Y,Z,In);
shading interp
colormap jet
axis tight
hold on
% plot3(0,0,0);

%
% view([az_tgt,el_tgt])

error('DONE');

%% calculate intensity

% renormalization
[Alpha,Beta] = meshgrid(dalpha,dbeta);
Gamma = sqrt(1-Alpha.^2 - Beta.^2);
g_real = imag(Gamma)==0;
L_tot = trapz(dbeta,trapz(dalpha,double(L),2))
L_real = trapz(dbeta,trapz(dalpha,double(L).*g_real,2));
K = L_tot/L_real
% U_tot = double(sum(abs(L(:)).^2));
% U_real = double(sum(abs(L(g_real)).^2));
% K = U_tot/U_real
L_p = K*L.*g_real; %*gamma_i*(illumination_params.lambda^2)/A_s
L_p_tot = trapz(dbeta,trapz(dalpha,double(L_p),2))

max_L = max(L_p(:))

I = L_p.*Gamma.*A_s;
I_tot = sum(I(:))
max_I = max(I(:))

I = I./0.2975416;


% renormalization
% U = abs(U).^2;
% L_tot = trapz(dbeta,trapz(dalpha,U,2));
% L_real = trapz(dbeta,trapz(dalpha,U.*g_real,2));
% K = L_tot/L_real
% % U_tot = double(sum(abs(L(:)).^2));
% % U_real = double(sum(abs(L(g_real)).^2));
% % K = U_tot/U_real
% A_s = physical_params.pitch*(physical_params.N*physical_params.M)*physical_params.ff;% area of diffracting aperature
% L_p = K*abs(U).^2.*g_real; %*gamma_i*(illumination_params.lambda^2)/A_s

% I_old = abs(U_old).^2;

%% plot pattern
dcs2ang = @(x) atand(tan(x));
% dcs2ang = @(x) deg2rad(atand(tan(x)));
% [Alpha,Beta] = meshgrid(dalpha,dbeta);

DI = figure('Name','Diffracted Irradiance');
% surf(Alpha,Beta,nthroot(I,1),'EdgeColor','none');
% view(2);
imagesc(dalpha,dbeta,nthroot(I,1));
% imagesc(dalpha,dbeta,nthroot(L_p,4));
% imagesc(dalpha,dbeta,log(L_p));
% imagesc(dcs2ang(dalpha),dcs2ang(dbeta),nthroot(L_p,4));
% imagesc(dcs2ang(dalpha),dcs2ang(dbeta),nthroot(L_p,4));
xlabel('\alpha')
ylabel('\beta')
% imagesc(dcs2ang(dalpha),dcs2ang(dbeta),I);
% surf(dcs2ang(Alpha),dcs2ang(Beta),nthroot(I,1),'EdgeColor','none');
set(gca,'YDir','normal');
% xlabel('Tilt Angle (\circ)')
% ylabel('Tip Angle (\circ)')
% zlabel('I')
daspect([1,1,1])
colormap hot
cb = colorbar;
% cb.Label.String = 'Diffracted Intensity (W/steradian projected area)';
hold on
plot(cos(0:pi/100:2*pi),sin(0:pi/100:2*pi))

data = struct('physical_params',pparams,...
              'illumination_params',illumination_params,...
              'bitmask',bitmask,...
              'view_params',view_params,...
              'alpha_axis',dalpha,...
              'beta_axis',dbeta,...
              'L',L);

ts = datetime('now','Format','HH-mm-ss_dd-MM-yyyy');
save([char(ts),config_name, '.mat'],'data');

threshold = 0.01;

binarized = I>(threshold*max(I(:)));
figure('Name','Binarized Locations'), imagesc(dalpha,dbeta,binarized),set(gca,'YDir','normal');
xlabel('\alpha'),ylabel('\beta'), colormap hot, colorbar, daspect([1 1 1])

stats = regionprops(binarized,'Centroid','BoundingBox','Area','Eccentricity','PixelList')
hold on

figure(DI);
for i = 1:length(stats)
    hold on
    bb = stats(i).BoundingBox;
    bbdcs = [interp1(dalpha,bb(1)),interp1(dbeta,bb(2)),abs(diff(interp1(dalpha,[bb(3)+bb(1),bb(1)]))),abs(diff(interp1(dbeta,[bb(4)+bb(2),bb(2)])))];
    rectangle('Position', bbdcs, 'EdgeColor','g','LineStyle',':')
    hold on
    scatter(stats(i).Centroid(1),stats(i).Centroid(2))
end


%% construct all alpha,beta needed here - we don't need to sample a uniform grid here, and can only sample where wethink there is light
alpha_sample_density = 1024/0.05;
beta_sample_density = 1024/0.05;
max_samples = 1024;

% Alphas = cell(1,length(stats));
% Betas  = cell(1,length(stats));

[a_min, a_max, b_min, b_max, s_samples] = deal(zeros(1,length(stats)));
[dalphas,dbetas,Ls] = deal({});



N = 731; % mirror rows - for now, enfore Odd
M = 731;% mirror cols - for now, enfore Odd
% N = 111; % mirror rows - for now, enfore Odd
% M = 111;% mirror cols - for now, enfore Odd
% N = 731; % mirror rows - for now, enfore Odd
% M = 731;% mirror cols - for now, enfore Odd
% % rotation params - not yet implemented
% phi_ipd = 0; % in plane rotation (degrees)
%% physical mirror params
PITCH = 13.68e-6; % mirror pitch (m)
FF = 0.92; % mirror fill factor
% R_m = @(lambda) 0.96; % mirror refelctivitiy function
R_m = 0.96*ones(N,M);
% H = 0e-6; % mirror height (m)
H = FF*PITCH*sind(12)*ones(N,M);
theta_dist = 'gaussian'; % 
theta_nom = 12; % nominal tilt angle 
theta_std = 0; % tilt angle stddev
%% substrate params
R_s = @(x) 0.6; % Substrate refelctivitiy function
%% Illumination Params
lambda = 4e-6; % wavelength (m)
theta_i = 24; % incidence elevation angle (deg)
% tehta_i = 5.956653; % optimal calculated from blaze condition
phi_i = 45; % incidence azimuth angle (deg)

%% Dervived params
A_s = N*M*PITCH;
% build mirror dist
switch theta_dist
    case 'gaussian'
        m_theta = (theta_std^2*(randn(N,M))+theta_nom);
    otherwise
        error('Unrecognized tilt angle distrubition')
end
%% bitmask patterning
% load('bitmask.mat')
% bitmask = data.bitmask;
% bitmask = ones(N,M);
% bitmask = -ones(N,M);
% bitmask = zeros(N,M);
% bitmask(2:2:end,:) = -1;
% bitmask(:,2:2:end) = -1;
bitmask = checkerboard(1,N,M)>0.5;
bitmask = bitmask(1:N,1:M); 
% % % bitmask = randi([0 1], N,M);
% bitmask = 2*bitmask - 1;
% bitmask = -1*bitmask;

r = 1;
% dc_alpha = 0.2;
% dc_beta = 0.2;
% coord = [r*dc_alpha, r*dc_beta, sqrt(r^2 - r^2*(dc_alpha^2 + dc_beta^2))];
% r = coord(3);
% bitmask = generate_zone_plate(lambda,theta_nom,N,M, coord);
% bitmask = bitmask';
% bitmask = -1*bitmask;

figure('Name','Bitmask'), imagesc(bitmask), colormap gray, caxis([-1 1]), colorbar,set(gca,'YDir','normal')
% pause
%% build param structs
% mirror array params
pp = struct(...
    'pitch',PITCH,...
    'ff',FF,...
    'R_m',R_m,...
    'R_s',R_s,...
    'N',N,...
    'M',M,...
    'h',H,...
    'm_theta',m_theta,...
    'ip_rot',0 ...
);

n_start = ((-(N-1)/2));
m_start = ((-(M-1)/2));
[mirror_n, mirror_m] = meshgrid((1:N)+n_start-1,(1:M)+m_start-1);
center_n = 0;
center_m = 0;
beam_waist_n = N/2; % in mirrors
beam_waist_m = M/2; % in mirrors

R = [cosd(phi_i) -sind(phi_i); sind(phi_i) cosd(phi_i)];
mirrors = R*[(mirror_n(:)-center_n)';(mirror_m(:)-center_m)'];
mirrors_n = reshape(mirrors(1,:),N,M);
mirrors_m = reshape(mirrors(2,:),N,M);

% E_m = exp(-pi*(((mirror_n-center_n)/(beam_waist_n)).^2 + ((mirror_m-center_m)/(beam_waist_m)).^2));
E_m = exp(-pi*(((mirrors_n)/(beam_waist_n)).^2 + (((mirrors_m)/(beam_waist_m)).^2).*(secd(theta_i)).^2));
figure('Name','Illumination Intensity'), imagesc(E_m), colorbar
% figure('Name','Illumination Intensity SURF'), surf(mirror_n,mirror_m,E_m,'EdgeColor','none'), colorbar
E_src = sum(E_m(:))* A_s;
I_src = E_src.^2 ;


% window_size = 20;
% waist = 0.5;
% [E_m,e_theta] = gaussian_beam(theta_i,phi_i, waist, 3*N*window_size);
% figure, imagesc(E_m), colorbar
% figure, surf(e_theta,e_theta,E_m,'EdgeColor','none'), colorbar
% E_m = conv((1/(window_size^2))*ones(window_size),E_m,'same')

% illumination params
ip = struct(...
    'lambda',lambda,...
    'theta',theta_i,...
    'phi',phi_i,...
    'center_n',center_n,...
    'center_m',center_m,...
    'beam_waist_n', beam_waist_n,...
    'beam_waist_m', beam_waist_m,... 
    'E_m',E_m...
);


for i = 1:length(stats)
    hold on
    bb = stats(i).BoundingBox;
    bbdcs = [interp1(dalpha,bb(1)),interp1(dbeta,bb(2)),abs(diff(interp1(dalpha,[bb(3)+bb(1),bb(1)]))),abs(diff(interp1(dbeta,[bb(4)+bb(2),bb(2)])))];
    corners = bbox2points(bbdcs);
    a_min(i) = min(corners(:,1));
    a_max(i) = max(corners(:,1));
    b_min(i) = min(corners(:,2));
    b_max(i) = max(corners(:,2));
    alpha_extent = a_max(i) - a_min(i);
    beta_extent = b_max(i) - b_min(i);
    box_samples = 2.^nextpow2(arrayfun(@(x) min(x,max_samples),[alpha_extent*alpha_sample_density,beta_extent*beta_sample_density]));
    s_samples(i) = max(box_samples);
    [Alphas{i},Betas{i}] = meshgrid(linspace(a_min(i),a_max(i),box_samples(1)),linspace(b_min(i),b_max(i),box_samples(2)));
%     centroid = stats(i).Centroid;
% end
    vp = struct(...
        'view_style','abs',... % 'abs','order','full'
        'r',r,...
        'alpha_min',a_min(i),...
        'alpha_max',a_max(i),...
        'beta_min',b_min(i),... 
        'beta_max',b_max(i),...
        'units','dcs',... % 'dcs','rad','deg'
        'samples',s_samples(i) ... prefer power of 2, max 1024 right now, need to subdivide grid in code for handling denser sampline since it is tied to cuda threadblock size
    );

    %% calculate field and intensity pattern
    fprintf([num2str(i) '/' num2str(length(stats)) ':  '])
    tic
    [dalphas{i},dbetas{i},Ls{i}] = dmd_diffracted_2(pp,ip,bitmask,vp,false);
    toc
end
%% plot 
figure('Name','Foveated Rendering')
for i = 1:length(stats)
    hold on
    surf(dalphas{i},dbetas{i},Ls{i},'EdgeColor','none')
    hold on
    bb = stats(i).BoundingBox;
    bbdcs = [interp1(dalpha,bb(1)),interp1(dbeta,bb(2)),abs(diff(interp1(dalpha,[bb(3)+bb(1),bb(1)]))),abs(diff(interp1(dbeta,[bb(4)+bb(2),bb(2)])))];
    rectangle('Position', bbdcs, 'EdgeColor','g','LineStyle',':')
end
grid off
% grid on
set(gca,'Color','k')
daspect([1 1 1])
colormap hot, 
colorbar,
view(2)

% % flattens everthing - would need to chunk out again...
% CAlpha = cell2mat(cellfun(@(x) reshape(x,1,[]), Alphas, 'UniformOutput', false));
% CBeta = cell2mat(cellfun(@(x) reshape(x,1,[]), Betas, 'UniformOutput', false));
% CLs = cell2mat(cellfun(@(x) reshape(x,1,[]), Ls, 'UniformOutput', false));
% CLs = reshape(CLs,length(CAlpha),length(CBeta));
% 
% figure, surf(CAlpha,CBeta,CLs,'EdgeColor','none')

%% old
% figure, 
% % surf(Alpha,Beta,nthroot(I,1),'EdgeColor','none');
% % view(2);
% imagesc(dalpha,dbeta,nthroot(I,5));
% xlabel('\alpha')
% ylabel('\beta')
% % imagesc(dcs2ang(dalpha),dcs2ang(dbeta),I);
% % surf(dcs2ang(Alpha),dcs2ang(Beta),nthroot(I,1),'EdgeColor','none');
% set(gca,'YDir','normal');
% % xlabel('Tilt Angle (\circ)')
% % ylabel('Tip Angle (\circ)')
% % zlabel('I')
% daspect([1,1,1])
% colormap hot
% cb = colorbar;
% % cb.Label.String = 'Diffracted Intensity (W/steradian projected area)';

% figure,
% % imagesc(dcs2ang(dalpha),dcs2ang(dbeta),I);
% % imagesc(dalpha,dbeta,I);
% % set(gca,'YDir','normal');
% surf(Alpha,Beta,nthroot(I_old,2),'EdgeColor','none');
% view(2);
% xlabel('\alpha')
% ylabel('\beta')
% % zlabel('I')
% daspect([1,1,1])
% % colormap hot
% colorbar

% err = abs(I2-I)./I*100;
% stats = datastats(err(:))
% figure, 
% surf(Alpha,Beta,err,'EdgeColor','none');
% view(2);
% xlabel('\alpha')
% ylabel('\beta')
% % zlabel('I')
% daspect([1,1,1])
% cb = colorbar;
% cb.Label.String = 'Error (%)';
% 
% figure, 
% histogram(err,50,'Normalization','probability')