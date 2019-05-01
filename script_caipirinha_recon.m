%%-----------------------------------------
%% SET RECONSTRUCTION PARAMETERS
%%-----------------------------------------

set(0,'DefaultFigureWindowStyle','docked')

img_size = [240,234,198];           % size of final reconstruction
Rtot = 9;                           % total acceleration factor
Ry = 3;                             % acceleration along Ry
del = 1;                            % CAIPIRINHA shift
num_chan = 20;                      % number of coil channels


%%-----------------------------------------
%% LOAD RAW DATA
%%-----------------------------------------

load data_caipi                     % load raw data for CAIPIRINHA
load sens_white                     % load coil sensitivity map (ESPIRiT)

imagesc3d2(rsos(data_caipi,4), [80,120,100], 11, [180,0,0], [0,0.3e-7], 0, 'Aliased CAIPIRINHA data at R=3x3')


%%-----------------------------------------
%% CAIPIRINHA: CALCULATE ALIASING CELL AND COLLAPSING INDICES
%%-----------------------------------------

im_size = img_size(2:3);

Rz = Rtot / Ry;

if del > Rz-1
    del = Rz-1;
end

mask_cell = zeros(Rtot);
mask_cell(1:Ry:end, 1:Rz:end) = 1;

for t = 1 : Ry : Rtot
    mask_cell(t,:) = circshift( mask_cell(t,:), [0, del * (t-1) / Ry] );
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


%%-----------------------------------------
%% CAIPIRINHA RECONSTRUCTION
%%-----------------------------------------


size_x = img_size(1);

i_ind = zeros(num_chan*Rtot*size_x, 1);
for t = 1:size_x
    ind = 1+(t-1)*num_chan : t*num_chan;
    ind = repmat(ind, [Rtot,1]);
    i_ind(1+(t-1)*Rtot*num_chan : t*Rtot*num_chan) = ind(:);
end

j_ind = zeros(num_chan*Rtot*size_x, 1);
for t = 1:size_x
    ind = 1 + (t-1)*Rtot : t*Rtot;
    ind = repmat(ind(:), [1,num_chan]);    
    j_ind(1+(t-1)*Rtot*num_chan : t*Rtot*num_chan) = ind(:);
end


sp_eye = speye(Rtot * size_x) * 1e-6;
g_factor = zeros(img_size);
    
% cov matrix identity, assume prewhitened sens maps
covmtx = eye(num_chan);
covmtx_h = sparse(kron(eye(size_x), covmtx));
covmtxinv_h = inv( covmtx_h );


img_caipi = zeros(img_size);
warning('off','all');

rec = permute(sens_white, [2,3,4,1]);

tic
for zi = 1:img_size(3)/Rz

    disp(num2str(zi))
    
    for cey = 1:img_size(2)/Ry
        % collapsing indices
        zi_ind = z_ind(cey,zi,:);
        cey_ind = y_ind(cey,zi,:);

        encoding = zeros(Rtot, num_chan, size_x);
        for l = 1:Rtot
            encoding(l,:,:) = rec(cey_ind(l),zi_ind(l),:,:);
        end

        data = permute(squeeze(data_caipi(:,cey,zi,:)), [2,1]);
                
        data = double(data(:));

        sparse_encoding_matrix = sparse(i_ind, j_ind, encoding(:), num_chan*size_x, Rtot*size_x, num_chan*Rtot*size_x);
                    
       
        res = sparse_encoding_matrix \ data;
        res = permute(reshape(res, [Rtot, size_x]), [2,1]);

        for l = 1:Rtot
            img_caipi(:,cey_ind(l),zi_ind(l)) = res(:,l);                 
        end            
                           
     end
end
toc  
 

imagesc3d2(img_caipi, [80,120,100], 21, [180,0,0], [0,0.7e-8], 0, 'CAIPIRINHA at R=3x3 with shift 1')





