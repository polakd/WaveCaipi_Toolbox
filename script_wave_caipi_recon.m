%%-----------------------------------------
%% SET RECONSTRUCTION PARAMETERS
%%-----------------------------------------

set(0,'DefaultFigureWindowStyle','docked')

img_size = [240,234,198];           % size of final reconstruction
Rtot = 9;                           % total acceleration factor
Ry = 3;                             % acceleration along Ry
del = 1;                            % CAIPIRINHA shift
num_chan = 20;                      % number of coil channels
OS = 6;                             % oversampling along the readout

%%-----------------------------------------
%% LOAD RAW DATA
%%-----------------------------------------

load data_wave                      % load raw data for Wave Caipi
load sens_white                     % load coil sensitivity map (ESPIRiT)
load psf                            % load PSF for Wave

mask_roi = ((sum(abs(sens_white),4)) ~= 0);
psf_opt = repmat(psfy,[1,1,img_size(3)]).*repmat(permute(psfz,[1,3,2]),[1,img_size(2),1]); 

imagesc3d2(rsos(data_wave,4), [720,120,100], 1, [180,0,0], [0,0.7e-8], 0, 'Aliased Wave data at R=3x3')


%%-----------------------------------------
%% WAVE CAIPI: CALCULATE ALIASING CELL AND COLLAPSING INDICES
%%-----------------------------------------

im_size = img_size(2:3);


Rz = Rtot / Ry;

if del > Rz-1
    del = Rz-1;
end

mask_cell = zeros(Rtot);
mask_cell(1:Ry:end, 1:Rz:end) = 1;

for t = 1 : Ry : Rtot
    mask_cell(t,:) = circshift( mask_cell(t,:), [1, del * (t-1) / Ry] );
end


alias_cell = fft2( mask_cell ) / Rtot;

[y_hat, z_hat] = find( alias_cell > 1e-2 );

[z, y] = meshgrid( 1:im_size(2), 1:im_size(1) );

y_ind = zeros([im_size, Rtot]);
z_ind = zeros([im_size, Rtot]);

for l = 1:Rtot
    y_ind(:,:,l) = y + (y_hat(l) - 1) * im_size(1) / Rtot;
    z_ind(:,:,l) = z + (z_hat(l) - 1) * im_size(2) / Rtot;
end

y_ind = y_ind .* (y_ind <= im_size(1)) + ( y_ind - im_size(1) ) .* ( y_ind > im_size(1) );
z_ind = z_ind .* (z_ind <= im_size(2)) + ( z_ind - im_size(2) ) .* ( z_ind > im_size(2) );


data_wave_fft = fftshift(fft(fftshift(data_wave,1),[],1),1) / sqrt(size(data_wave,1));



%%-----------------------------------------
%% Wave CAIPI RECONSTRUCTION
%%-----------------------------------------

img_wave = zeros(img_size);

lsqr_tol = 1e-3;
lsqr_iter = 300;
    
warning('off','all');

param = [];
param.psf_len = img_size(1)*OS;
param.img_len = img_size(1);
param.pad_size = img_size(1)*(OS-1)/2;
param.num_chan = num_chan;
param.Rtot = Rtot;


tic
for zi = 1:img_size(3)/Rz
%for zi = img_size(3)/Rz/2

    disp(num2str(zi))
    
    for cey = 1:img_size(2)/Ry
        % collapsing indices

        zi_ind = z_ind(cey,zi,:);
        cey_ind = y_ind(cey, zi,:);
        
        
        psfs = zeros([img_size(1)*OS,Rtot,num_chan]);
        
        mask_use = zeros([img_size(1),Rtot]);
        rcv = zeros([img_size(1),Rtot,num_chan]);

        for l = 1:Rtot
            mask_use(:,l) = mask_roi (:, cey_ind(l),zi_ind(l));
            psfs(:,l,:) = repmat( psf_opt(:,cey_ind(l),zi_ind(l)), [1,1,num_chan] );
            rcv(:,l,:) = sens_white(:,cey_ind(l),zi_ind(l),:);
        end

        if sum(sum(sum(mask_use,1),2),3) > 0
        
            rhs = squeeze( data_wave_fft(:,cey,zi,:) );

            param.psfs = psfs;
            param.rcv = rcv;

            [Res,~] = lsqr(@apply_wave_caipi, rhs(:), lsqr_tol, lsqr_iter, [], [], [], param);        

            Res = reshape(Res, [param.img_len, Rtot]);

            for l = 1:Rtot
                img_wave(:,cey_ind(l),zi_ind(l)) = Res(:,l);                 
            end
            
        end
        

     end
end
toc  


imagesc3d2(img_wave, [80,120,100], 2, [180,0,0], [0,0.7e-8], 0, 'Wave-Caipi at R=3x3 with shift 1')


