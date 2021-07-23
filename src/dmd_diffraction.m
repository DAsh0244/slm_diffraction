function [dalpha,dbeta,U] = dmd_diffraction(physical_params,illumination_params,bitmask,view_params,use_double)
%dmd_diffraction Calculates the Diffracted field from a DMD
% Usage:
%   
%   [alpha,beta,U] = dmd_diffraction(physical_params,illumination_params,bitmask,viewing_params);
% 
% Input Parameter Description:
%   physical_params is a struct containing attributes:
%       pitch: the physical mirror pitch (m)
%       ff: array fill factor (1.0 is 100% fill factor)
%       R_m: mirror reflectance function taking illumination wavelength as input arguement
%       R_s: substrate reflectance function taking illumination wavelength as input arguement
%       N: Number of rows in array (integer value)
%       M: Number of columns in array (integer value)
%       h: mirror height (m) - This is a pure piston term globally applied
%       m_theta: mirror tilt distrubition array. In normal operation this
%                is all a single value of the nominal tilt angle possible for the mirrors, 
%                and takes the form of either (1) an (NxM) array of tilt angles for each mirror, 
%                or (2) a (1x(N*M)) or ((N*M)x1) vector of tilt angles for each mirror.
%                The values here are then modulated by the bitmask parameter.
%                For further detail refer t oteh bitmask parameter descroption. 
%       ip_rot: in-plane rotation angle (degrees) - NOT YET IMPLEMENTED
%
%   illumination_params is a struct containing attributes:
%         lambda: wavelength (m)
%         theta: incident elevation angle from substrate normal (degrees). Ranges from [0,90] deg
%         phi:  incident azimuthal angle (degrees). Ranges from [0,359] deg
%         E_m:  incident field amplitude at every mirror
%
%   bitmask is an array that describes mirror states. 
%       The array takes the form of an (NxM) array or (1x(N*M)) or ((N*M)x1) vector of values to describe mirror states. 
%       Typically is made of 1/-1 elements to represent +/- tilt angle operation.
%       Values in the range [-1,1] are acceptable, although integer numbers are the most realistic to use.
%       0 values would imply flat state. The values here are used to elementwise multiply
%       the tilt angles provided in physical_params.m_theta and thuse must be the same shape arrays. 
%
%   view_params is a struct containing attributes:
%       view_style:
%       r:
%       alpha_min: 
%       alpha_max: 
%       beta_min: 
%       beta_max:
%       order_x:
%       order_y:
%       spacing_param:
%       samples: 
%
% Return parameter Description:
%  

%% DEV zone - tweaking parameters here for performance reasons 
% (make the final 'use_double' a struct of parameters for more variable
% tuning? eg: <CPU|GPU>, 
%             <double|single>,
%             <sinc|LUT>, specify the sinc implemetation externally and pass a function handle?
%             <memory_margin - Althought should be calucaltable with some more thought>
%             <matlab|C> calculation?

  mem_margin = 20; % how much to pad memory requirements to prevent OOM errors
  sinc_implementation = @sinc; % implementation of a vectorized sinc function taking complex floating values as input. 

  % setup numerical types used for calcuation and memory requirements
  if use_double
    numerical_type = @double;
    bytes_per_entry = 8;
  else
    numerical_type = @single;
    bytes_per_entry = 4;
  end

  % ensure gpu is present
  assert(gpuDeviceCount ~= 0,'NO GPU FOUND!')

  % grab GPU handle
  gpu = gpuDevice();

%% dervived params
% mirror params
% gap = (1 - physical_params.ff) * physical_params.pitch / 2; % unused
  w = physical_params.ff*physical_params.pitch;
  num_mirrors = physical_params.N*physical_params.M;
  if isa(physical_params.R_m,'function_handle')
    R_m = physical_params.R_m(illumination_params.lambda,physical_params.N,physical_params.M);
  else % assume reflectivities are given in a physical_params.N physical_params.M matrix or equivalanet row/col vector;
    R_m = physical_params.R_m;
  end

% normalized params
  delta = physical_params.pitch ; % normalized pitch
  w_hat = w / illumination_params.lambda; % normalized mirror width
  h = physical_params.h;%  / illumination_params.lambda; % normalized mirror height
  tantheta = tand(physical_params.m_theta.*bitmask); % precalculate tan(mirror_angles)

% incidence angle params - modified by PR
  alpha_i = illumination_params.alpha_i; 
  beta_i  = illumination_params.beta_i; 
  gamma_i = sqrt( 1 - (alpha_i.^2+beta_i.^2) );

% viewing space params
%   r_hat = view_params.r/illumination_params.lambda; % normalized viewing radius
  r = view_params.r;
  lambda = illumination_params.lambda;
  alpha_target = -alpha_i+(view_params.order_x*illumination_params.lambda/physical_params.pitch);
  beta_target = -beta_i+(view_params.order_y*illumination_params.lambda/physical_params.pitch);
  if (length(view_params.samples) >1)
    samples_x = view_params.samples(1);
    samples_y = view_params.samples(2);
  else
    samples_x = view_params.samples;
    samples_y = view_params.samples;
  end
  
  switch view_params.view_style
    case 'abs' % absolute coordinates
      dalpha = linspace(view_params.alpha_min,view_params.alpha_max,samples_x);
      dbeta  = linspace(view_params.beta_min,view_params.beta_max,samples_y);
    case 'order' % target a region around a specific diffracted order
      alpha_target = -alpha_i+(view_params.order_x*illumination_params.lambda/physical_params.pitch);
      beta_target = -beta_i+(view_params.order_y*illumination_params.lambda/physical_params.pitch);
      dalpha = linspace(alpha_target-view_params.spacing_param,alpha_target+view_params.spacing_param,samples_x);
      dbeta  = linspace(beta_target-view_params.spacing_param,beta_target+view_params.spacing_param,samples_y);
    case 'full' % sample the whole sphere
      dalpha = linspace(-1,1,max(view_params.samples));
      dbeta = linspace(-1,1,max(view_params.samples));
    case 'evanescent'
      dalpha = [linspace(-2,-1,max(view_params.samples)),linspace(1,2,max(view_params.samples))];
      dbeta = linspace(-2,2,max(view_params.samples));
    otherwise
      error('Unrecognized view_style')
  end
  [Alpha,Beta]= meshgrid(dalpha,dbeta);
  % Modified by PR - capitalized Gamma
  Gamma = sqrt(1-Alpha.^2 - Beta.^2);
  g_real = imag(Gamma)==0;

%% Calculate field - modified by PR
% This is done on a chunk basis, where the data cube is defined as having
% axes alpha,beta,mirror_idx. From each chunk, it is compressed along
% the mirror_idx via summation. This process is repeated until all mirrors
% have been calculated.

% generate n,m values coresponding to mirrors placements
% couldprpobably get away with reshaping a meshgrid instead of this...
  parity = mod([physical_params.N,physical_params.M],2);
  if parity(1)==1 % Odd array
      n_start = ((-(physical_params.N-1)/2));
      n_list = (1:physical_params.N)+n_start-1;
  else % even array
      % build array by building positive half, then mirror.
      n_list = (1:physical_params.N/2) - 1/2;
      n_list = [-flip(n_list),n_list];
  end
  if parity(2) == 1 % Odd array
      m_start = ((-(physical_params.M-1)/2));
      m_list = (1:physical_params.M)+m_start-1;
  else % even array
      m_list = (1:physical_params.M/2) - 1/2;
      m_list = [-flip(m_list),m_list];
  end
  
  
% n_vals = gpuArray(reshape(numerical_type(idivide(int32(1:num_mirrors),int32(physical_params.N),'ceil')+n_start-1),[1 1 num_mirrors]));
% m_vals = gpuArray(reshape(numerical_type((mod((1:num_mirrors)-1,physical_params.M)+1)+m_start-1),[1 1 num_mirrors]));
  [n,m] = meshgrid(n_list,m_list);

  n_vals = gpuArray(reshape(numerical_type(n),[1 1 num_mirrors]));
  m_vals = gpuArray(reshape(numerical_type(m),[1 1 num_mirrors]));
  %
  mirror_tilts    = gpuArray(reshape(numerical_type(tantheta),[1 1 num_mirrors]));
  % UNUSED FOR NOW
  mirror_heights  = gpuArray(reshape(numerical_type(h),[1 1 num_mirrors])); % - NOT IMPLEMENTED
  %
  mirror_reflectivities = gpuArray(reshape(numerical_type(R_m),[1 1 num_mirrors]));
  mirror_intensities    = gpuArray(reshape(numerical_type(illumination_params.E_m),[1 1 num_mirrors]));
  mirror_intensities = mirror_intensities .* mirror_reflectivities;
  clear mirror_reflectivities

  Alpha = gpuArray(Alpha);
  Beta  = gpuArray(Beta);
  Gamma  = gpuArray(Gamma);

%% chunking
  BLOCK_SIZE = gpu.MaxThreadBlockSize;
  % Modified by PR
%   BLOCK_SIZE = [max(BLOCK_SIZE(2),samples_y),min(BLOCK_SIZE(1),samples_x),min(BLOCK_SIZE(3)/2,num_mirrors)];
  BLOCK_SIZE = [samples_y,samples_x,min([BLOCK_SIZE(3)/4,num_mirrors])];
  NUM_BLOCKS = ceil(num_mirrors/BLOCK_SIZE(3));
  BLOCK_MEM = 2*bytes_per_entry*prod(BLOCK_SIZE);
  BLOCKS_PER_CHUNK = max(1,floor(gpu.AvailableMemory/(mem_margin*BLOCK_MEM))); % mem_margin is a fudge factor
  BLOCKS_PER_CHUNK = 1;
  NUM_CHUNKS = ceil(NUM_BLOCKS/BLOCKS_PER_CHUNK);
  CHUNK_SIZE = BLOCK_SIZE.*[1 1 BLOCKS_PER_CHUNK];

%   matObj = matfile('phasors.mat','Writable',true);
%   matObj.phasors = zeros(BLOCK_SIZE(1), BLOCK_SIZE(2), num_mirrors);

  U = zeros(BLOCK_SIZE(1),BLOCK_SIZE(2),1, functions(numerical_type).function,'gpuArray');
%     progressbar('Chunks')
  textprogressbar('Calculating Chunks: ');
  for chunk = 0:NUM_CHUNKS-1
    % identify mirror indices
    idxs = int32(1+((chunk*CHUNK_SIZE(3)):min((chunk*CHUNK_SIZE(3))+CHUNK_SIZE(3)-1,num_mirrors-1))); % generate indicies for which slices of mirror values to use.

    % calculate ray lengths
    lx = bsxfun(@minus,r*Alpha,double(n_vals(idxs))*delta);
    ly = bsxfun(@minus,r*Beta,double(m_vals(idxs))*delta);
%     lz = bsxfun(@minus,r*Gamma,mirror_heights(idxs));
    lz = repmat(r*Gamma,[1 1 length(idxs)]);
    L_mn = sqrt(lx.^2 + ly.^2 +lz.^2);
    lx = lx ./L_mn;
    ly = ly ./L_mn;
    lz = lz ./L_mn;
        
    % alpha,beta,gamma,eff are modulated by the l vector
    alpha_eff = gpuArray(numerical_type(alpha_i + lx ));
    beta_eff  = gpuArray(numerical_type(beta_i  + ly ));
    gamma_eff = gpuArray(numerical_type(gamma_i + lz ));
    
    clear lx ly lz
    
%     prop_scalar = double(w^2*lambda*r/1i); % small number in physical coordinates, large when scaled
%     prop_term = double(bsxfun(@rdivide,Gamma,L_mn.^2)); % large number in both coords, really large in scaled
%     prop_amplitude = numerical_type(prop_scalar.*prop_term);
    prop_amplitude = numerical_type(...
        (w^2*lambda*r/1i).*bsxfun(@rdivide,Gamma,L_mn.^2)...
       );
    
%    L_mn = numerical_type(L_mn )
   
    % compute the tilt terms for computing the envelopes 
    tilt_term = bsxfun(@times, gamma_eff, ((mirror_tilts(idxs)./sqrt(2)))); 
    
    % unity phasors
%     spherical_wave = exp(1i*2*pi/lambda*L_mn); 
%     linear_phase = exp(-1i*2*pi/lambda*delta*(bsxfun(@times,alpha_i,n_vals(idxs))+bsxfun(@times,beta_i,m_vals(idxs))));
%     height_term = exp(-1i*2*pi/lambda*bsxfun(@times,gamma_eff,mirror_heights(idxs)));
%     
%     unity_phasors = linear_phase.*height_term.*spherical_wave;
  
    unity_phasors = ...
        exp(1i*2*pi/lambda*L_mn).*...
        exp(-1i*2*pi/lambda*delta*(bsxfun(@times,alpha_i,n_vals(idxs))+bsxfun(@times,beta_i,m_vals(idxs)))).*...
        exp(-1i*2*pi/lambda*bsxfun(@times,gamma_eff,mirror_heights(idxs)));
    
    clear prop_scalar prop_term 
    % Envelope cobtribution
    % ---------------------
    envelope = ...
        sinc_implementation(w_hat*(alpha_eff-tilt_term)).*...
        sinc_implementation(w_hat*(beta_eff-tilt_term));
    
    scaled_field = unity_phasors.*bsxfun(@times, mirror_intensities(idxs), envelope);
    % Diffracted Field
%     cube = bsxfun(@times, mirror_intensities(idxs), prop_amplitude.*scaled_field); 
    cube = prop_amplitude.*scaled_field;
%       matObj.phasors(:,:,idxs)=cube;

    % Add up mirror contributions
    U = U + sum(cube,3); % sum along mirror axis
%     U_tmp= 0;
%     for adx = 1:length(idxs)
%         n_vals(adx)
%         m_vals(adx)
%         figure(10), imagesc(abs(prop_term(:,:,adx)).^2), colorbar
%         temp = squeeze(exp(1i*2*pi/lambda*bsxfun(@times,mirror_heights(adx),gamma_eff)));
%         figure(11), imagesc(abs(envelope(:,:,adx)).^2), colorbar
%         figure(11), imagesc(abs(U).^2), colorbar
%         pause(0.05)
% %         figure(14), imagesc(abs(alpha_eff(:,:,adx)).^2), colorbar
% %         figure(15), imagesc(abs(beta_eff(:,:,adx)).^2), colorbar
%     end
%     alpha_i
%     beta_i
%     alpha_target
%     beta_target        
%     figure(12), imagesc(dalpha,dbeta, abs(cube(:,:,adx)).^2), colorbar
%         set(gca,'YDir','normal');
%         hold on
%         scatter(alpha_target,beta_target,'r')
%         
% %         hold off
%         U_tmp = U_tmp + cube(:,:,adx);
%         figure(11)
%         imagesc(dalpha,dbeta, abs(U_tmp).^2), colorbar
%         set(gca,'YDir','normal');
%         title(sprintf('m=%d,n=%d',m_vals(adx),n_vals(adx)))

%         hold on
%         scatter(alpha_target,beta_target,'r')
%         hold off
%         pause(0.05)
%     end
%       progressbar(chunk/(NUM_CHUNKS-1));
    textprogressbar(max([100*chunk/(NUM_CHUNKS-1),1]))
%       wait(gpu) % ensure U has been properly updated - not needed 
  end
  U = squeeze(gather(U));
  textprogressbar(' done');

%   matObj.dalpha=dalpha;
%   matObj.dbeta=dbeta;
%   matObj.bitmask = bitmask;

end
