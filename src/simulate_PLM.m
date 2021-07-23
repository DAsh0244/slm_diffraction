clear all
% close all
clc
config_name = 'PLM';% 'optimization';

DIRECTORY = '..\data\';%'..\optimization\';



%% Mirror array params
N = 16*20/2;% mirror rows - for now, enfore Odd
M = 16*20/2;% mirror cols - for now, enfore Odd
NUM_BITS = 4;
NUM_LEVELS = 2^NUM_BITS;
levels = [0,0.01,0.0205,0.0422,0.0560,0.0727,0.1131,...
    0.1734,0.3426,0.3707,0.4228,0.4916,0.5994,...
    0.6671,0.7970,0.9375,1];
%% physical mirror params
PITCH = 10.8e-6; % mirror pitch (m)
FF = 10.5/10.8; %0.945; % mirror fill factor
% R_m = @(lambda) 0.96; % mirror refelctivitiy function
R_m = 0.92*ones(N,M);
w = FF*PITCH;
mirror_height = 0.532e-6/2; %0.243e-6;   %3*w/sqrt(2)*sind(12);

[X,Y] = meshgrid(0:M-1,0:N-1);


heightmap = round(15/(NUM_LEVELS-1)*mod(flip(X),NUM_LEVELS));
% heightmap = (cos(2*pi/NUM_LEVELS*(X-(N-1)/2))>0);
% NUM_LEVELS*(mod(X,NUM_LEVELS)<NUM_LEVELS/2);
% heightmap = zeros(N,M);
% heightmap(mod(0:N-1,2)==0,:) = 1;
% heightmap(:,mod(0:M-1,2)==0) = NUM_LEVELS;
% heightmap(:,1:2:end) = NUM_LEVELS;
figure('Name','Quantized Heightmap'),
imshow(heightmap,'InitialMagnification','fit');
colorbar; colormap gray;
caxis([0 NUM_LEVELS-1]);
set(gca,'YDir','normal');
title('Bit Mask DMD');

% bitmap =  gen_bit_map_p47(heightmap+1,'blazed_bitmap.bmp',false);
% imwrite(uint8(heightmap*255/15),'blazed_heightmapX.bmp')
% imwrite(bitmap, 'blazed_bitmapX.bmp')
% fit to real heights
% heightmap = levels(heightmap+1);
heightmap = heightmap./(NUM_LEVELS);
H = mirror_height .* heightmap; % loaded from MAT file

theta_dist = 'gaussian'; %
theta_nom = 0; % nominal tilt angle
theta_std = 0; % tilt angle stddev

%% substrate params
R_s = @(x) 0.6; % Substrate refelctivitiy function

%% Dervived params
A_s = N*M*PITCH*PITCH;
% build mirror dist
switch theta_dist
    case 'gaussian'
        m_theta = (theta_std^2*(randn(N,M))+theta_nom);
    otherwise
        error('Unrecognized tilt angle distrubition')
end

az_in = 0;
% az_in = -90;
steps = 30;
efficiencies = zeros(1,steps);
el_ins = 90-linspace(0,85,steps);
for index = 1:steps
    el_in = el_ins(index);
    % az_in = 45;       % x to y axis
    % el_in = (90-24);  % y to z axis
    uv_in = azel2uv([az_in;el_in]);
    % % dervived params
    beta_i = uv_in(1);
    gamma_i = uv_in(2)
    alpha_i = sqrt(1-(beta_i^2+gamma_i^2))
    
    %% Illumination Params
    lambda = 632e-9; % 0.532e-6; % wavelength (m)
    u = beta_i; % u = beta
    v = sqrt(1-(alpha_i.^2+beta_i.^2));     % v = gamma
    azel = uv2azel([u;v]);
    az_i = azel(1);  % x to y axis
    el_i = azel(2);  % y to z axis
    
    strBaseFileName = sprintf('PML_BIN_PHASE_%f_%f',az_i,el_i);
    
    fprintf('Incident Illumination - Azimuth,Elevation: %2.3f,%2.3f\n',az_i,el_i);
    % u = 0; %data.beta_tgt;                                % u = beta
    % v = 1;%sqrt(1-(data.alpha_tgt.^2+data.beta_tgt.^2)); % v = gamma
    % azel = uv2azel([u;v]);
    % az_tgt = azel(1);  % x to y axis
    % el_tgt = azel(2);  % y to z axis
    % fprintf('Target - Azimuth,Elevation: %2.3f,%2.3f\n',az_tgt,el_tgt);
    
    bitmask = 0*heightmap; %data.bitmask(1:M,1:N);
    
    %% build param structs
    % mirror array params
    pparams = struct(...
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
    [mirror_n, mirror_m] = meshgrid((1:M)+m_start-1,(1:N)+n_start-1);
    center_n = 0;
    center_m = 0;
    beam_waist_n = M/3; %365.4971/2;%N/2; % in mirrors
    beam_waist_m = M/3; %365.4971/2;%M/2; % in mirrors
    
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
    
    E_m = ones(size(E_m));
    
    % figure('Name','Illumination Intensity'),
    % imshow(E_m,'InitialMagnification','fit' );
    % set(gcf,'Units','Normalized');
    % colorbar
    % figPos = get(gcf,'OuterPosition');
    
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
        'alpha_i',alpha_i,...data.alpha_i,...
        'beta_i',beta_i,...data.beta_i,...
        'center_n',center_n,...
        'center_m',center_m,...
        'beam_waist_n', beam_waist_n,...
        'beam_waist_m', beam_waist_m,...
        'E_m',E_m...
        );
    
    %% View parameters
    FN = (M*PITCH)^2 / lambda
    r = FN*10; % data.R_0;       % Distance @ which to create spot
    % r = 2.5/100;       % Distance @ which to create spot
    
    % view_params = struct(...
    %     'view_style','order',... % 'abs','order','full'
    %     'r',r,...
    %     'order_x',1,...
    %     'order_y',0,...
    %     'spacing_param',0.05,...
    %     'units','dcs',... % 'dcs','rad','deg'
    %     'samples',1024 ... prefer power of 2, max 1024 right now, need to subdivide grid in code for handling denser sampline since it is tied to cuda threadblock size
    % );
    
    
    %

    % view_params = struct(...
    %     'view_style','full',... % 'abs','order','full'
    %     'r',r,...
    %     'order_x',0,...data.frac_order_p,...
    %     'order_y',0,...data.frac_order_q,...
    %     'spacing_param',0.025,...125
    %     'units','dcs',... % 'dcs','rad','deg'
    %     'samples',512 ... prefer power of 2, max 1024 right now, need to subdivide grid in code for handling denser sampline since it is tied to cuda threadblock size
    % );
    
    %% refernce img
%     if abs(alpha_i) > 0.45 
        view_params = struct(...
        'view_style','abs',... % 'abs','order','full'
        'r',r,...
        'beta_min',-0.05,...
        'beta_max',0.05,...
        'alpha_min',-1.7,...1,...25,...
        'alpha_max',1,...1,...25,...
        'order_x',-1/NUM_LEVELS,...
        'order_y',0,...
        'spacing_param',2e-3/2,...0.05/4,...
        'units','dcs',... % 'dcs','rad','deg'
        'samples',[4096*4,16] ... prefer power of 2, max 1024 right now, need to subdivide grid in code for handling denser sampline since it is tied to cuda threadblock size
        );
    
    tic
    [dalpharef,dbetaref,Uref] = ...
        dmd_diffraction(pparams,illumination_params,bitmask,view_params,false);
    toc
    Iref = abs(Uref).^2;
%     figure, imagesc(dalpharef,dbetaref,log10(Iref./max(Iref(:))))
%     colorbar
    [Alpha,Beta] = meshgrid(dalpharef,dbetaref);
    Gamma = sqrt(1-Alpha.^2 - Beta.^2);
    g_real = imag(Gamma)==0;
    
    
    U_tot = double(sum(Iref(:)));
    U_real = double(sum(Iref(g_real)));
    K = U_tot/U_real

    %% calculate field and intensity pattern
        view_params = struct(...
        'view_style','order',... % 'abs','order','full'
        'r',r,...
        'beta_min',-0.0002,...
        'beta_max',0.0002,...
        'alpha_min',-0.02,...1,...25,...
        'alpha_max',0.02,...1,...25,...
        'order_x',-1/NUM_LEVELS,...
        'order_y',0,...
        'spacing_param',2e-3/2,...0.05/4,...
        'units','dcs',... % 'dcs','rad','deg'
        'samples',[512*2,64] ... prefer power of 2, max 1024 right now, need to subdivide grid in code for handling denser sampline since it is tied to cuda threadblock size
        );
    
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
%     L_tot = trapz(dbeta,trapz(dalpha,double(L),2));
%     L_real = trapz(dbeta,trapz(dalpha,double(L).*g_real,2));
%     K = L_tot/L_real;
%     U_tot = double(sum(abs(Iref(:)).^2));
%     U_real = double(sum(abs(L(g_real)).^2));
%     K = U_tot/U_real
    % L_p = K*L.*g_real; %*gamma_i*(illumination_params.lambda^2)/A_s
    % L_p_tot = trapz(dbeta,trapz(dalpha,double(L_p),2));
    I = K*I;
    % max_L = max(L_p(:));
    
    % I_L = L_p.*Gamma.*A_s;
    
    % ISS = 2 * view_params.spacing_param ./ view_params.samples .* view_params.r;
    [nR,nC] = size(I);
    ISSX = (view_params.alpha_max - view_params.alpha_min  )./ view_params.samples(1) .* view_params.r;
    ISSY = (view_params.beta_max - view_params.beta_min  )./ view_params.samples(2) .* view_params.r;
    
    % ISSX = (view_params.spacing_param )./ view_params.samples(1) .* view_params.r;
    % ISSY = (view_params.spacing_param )./ view_params.samples(2) .* view_params.r;
    % ISSX = ISS;
    % ISSY = ISS;
    xg = -(nC-1)/2:(nC-1)/2;  xg = xg * ISSX;
    yg = -(nR-1)/2:(nR-1)/2;  yg = yg * ISSY;
    
    sI = g_real.*I./max(I(:));
    [m,max_index] = max(sI(:));
    [r,c] = ind2sub(size(sI),max_index);
    % figure, imagesc(dalpha,dbeta,sI), colormap hot, colorbar
    % figure, imagesc(xg/1e-3,yg/1e-3,nthroot(sI,1)), colormap hot, colorbar
    % daspect([1 1 1])
    % set(gca,'YDir','normal');
    
    % figure, plot(xg/1e-3, sI(r,:))
    % set(gca,'Yscale','log');
    
    % calculate order locations
    avg_sI = sum(sI);
    energy = sum(avg_sI);
    PERIOD = PITCH*NUM_LEVELS;
    fprintf('order at: %f',-alpha_i +(-1*lambda/PERIOD))
    % figure,
    % ax1 = subplot(211);
    % plot(dalpha, avg_sI./max(avg_sI)); %sI(r,:)) %)
    % orders_target = -5:5;
    % PERIOD = PITCH*2*8;
    % xlabel('Diffracted Angle (degrees)')
    % ylabel('Diffraction Effficiency (%)')
    %
    % for order=orders_target
    %     hold on;
    %     if order == 0
    %         xline(-alpha_i +(order*lambda/PERIOD),'r')
    %     else
    %         xline(-alpha_i +(order*lambda/PERIOD))
    %     end
    % end
    % ax2 = subplot(212);
    % % figure,
    % plot(dalpha, cumsum(100*avg_sI./energy)); %sI(r,:)) %)
    % linkaxes([ax1, ax2],'x')
    %
    % % figure,
    % % plot(dalpha, sI(r,:)); %avg_sI./max(avg_sI)); %sI(r,:)) %)
    %
    % for order=orders_target
    %     hold on;
    %     if order == 0
    %         xline(-alpha_i +(order*lambda/PERIOD),'r')
    %     else
    %         xline(-alpha_i +(order*lambda/PERIOD))
    %     end
    % end
    % beta_target = -beta_i+(view_params.order_y*illumination_params.lambda/physical_params.pitch);
    
    efficiencies(index) = sum(double(I(:)))/9.693607733986699e-25*100
    
%     figure, bar(1:NUM_LEVELS,2./lambda.*2.*H(1,(1:int32((NUM_LEVELS)))))
%     % figure, stairs(heightmap(1,1:int32(1*(NUM_LEVELS))))
%     title('height profile')
%     ylabel('Phase shift (x\pi)')
%     xlim([1 NUM_LEVELS])
    % ylim([0 2])
    % figure, plot(dbeta, sI(:,c))
    
    MAT_File = [DIRECTORY '\' strBaseFileName '.mat'];
    
    sim_res = struct('physical_params',pparams,...
        'illumination_params',illumination_params,...
        'bitmask',bitmask,...
        'heightmap',heightmap,...
        'view_params',view_params,...
        'alpha_axis',dalpha,...
        'beta_axis',dbeta,...
        'I_U',I,...
        'U',U);
    
    save(MAT_File,'sim_res','sim_res');
    
    
    % close all
end
% 2.509602445600622e-25*100 %1.274796855914139e-25*100
    figure, bar(1:NUM_LEVELS,2./lambda.*2.*H(1,(1:int32((NUM_LEVELS)))))
    % figure, stairs(heightmap(1,1:int32(1*(NUM_LEVELS))))
    title('height profile')
    ylabel('Phase shift (x\pi)')
    xlim([1 NUM_LEVELS])
    
figure, plot(90-el_ins,efficiencies)
aois = 90-el_ins;
figure, plot([-flip(aois) aois] ,[flip(efficiencies),efficiencies])
xlabel('Angle of Incidence (degreees)')
ylabel('Diffraction Efficiency (%)')
grid on
xlim([-90 90])