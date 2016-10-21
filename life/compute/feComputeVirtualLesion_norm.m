function [rmse_wVL, rmse_woVL, nFas_tract, nFas_PN, nVoxels] = ...
          feComputeVirtualLesion_norm(fe, ind_tract, varargin)
%
% This function compute the root-mean-squred error in the voxels defined by
% the anatomy of a white matter tract of interest (a set of fascicles
% identified by ind_tract).
%
% The RMSE is reutrned with (rmse_wVL) and without (rmse_woVL) the tract of
% interest.
%
% – INPUTS –
%   fe         : fe structure
%   ind_tract  : indices to fibers in the tract to be virtually lesioned
%
% – OUTPUTS –
%   rmse_wVL   : Root-mean-squared error of the lesioned model,
%                path-neighborhood fascicles only, no tract of interest.
%   rmse_woVL  : Root-mean-squared error of the full model with tract of
%                interest and path-neighborhood.
%   nFas_tract : Total number of fascicles in the tract of interest.
%   nFas_PN    : Total number of fascicles in the path-neighborhood of the
%                tract.
%   nVoxels    : Total number of voxels in the tract.
%
% 
%  Copyright (2016),   Franco Pestilli       | Cesar F. Caiafa
%  Indiana University, frakkopesto@gmail.com | ccaiafa@gmail.com

% Find the subset of voxels intersected by the tract of interest 
% (actually by all fascicles in the tract of interest).
%
% To do so, we will find all the entries in the third mode of the 
% tensor touched by the fascicles of the tract of interest.
[inds, ~] = find(fe.life.M.Phi(:,:,ind_tract));

% Then we will reduce the tensor to by eliminating the zero entries in the
% tensor Phi.
%
% The inds is a list of (i,j,k) positions of the nnz entries in Phi. 
% Because we are interested in locating voxels (mode 2) we look at 
% the second column of inds, j.
%
% The result is the set of indices of the voxels of the tract of interest.
voxel_ind = unique(inds(:,2));

% Now we find all fascicles crossing the voxels of the tract of interst.
% These contains both tract and path-neightborrhood fascicles.
%
% To do so we find the subtensor corresponding to the tract voxels.
inds_pn    = feGet(fe,'pathneighborhood',ind_tract);
nVoxels    = length(voxel_ind);
nTheta     = feGet(fe,'n bvals');

% Find the number of fascicles in the tract of interest and those in the
% path-neighborhood.
nFas_tract = length(ind_tract);
nFas_PN = length(inds_pn);

% COlelct the fascicles weights for the model with (w) and without (wo)
% virtual lesion (VL):
wo_VL = feGet(fe, 'fiber weights');
w_VL  = wo_VL;
w_VL(ind_tract) = 0;

% Extract the measured diffusion signal (anisotropic, or demeaned).
measured  = feGet(fe,'dsigdemeaned by voxel');

% Compute rmse restricted to Path neighborhood voxels without Virtual Lesion
S0       = fe.life.diffusion_S0_img(voxel_ind); % this shoudl change to an feGet call.
measured = measured(:,voxel_ind);

% Restrict tensor model to the voxels of the tract of interest.
M              =  fe.life.M; % This should an feGet call.
M.Phi          =  M.Phi(:,voxel_ind,:);

% Compute the RMSE of full, unlesioned, model.
predicted_woVL =  reshape(M_times_w(M.Phi.subs(:,1), ...
                                    M.Phi.subs(:,2), ...
                                    M.Phi.subs(:,3), ...
                                    M.Phi.vals, ...
                                    M.DictSig, ...
                                    wo_VL, ...
                                    nTheta, ...
                                    nVoxels), ...
                                    size(measured));

rmse_woVL = sqrt(mean((measured - predicted_woVL).^2,1));

% Compute the RMSE after removing the tract of interest fascicles 
% from the life model.
predicted_VL    = reshape(M_times_w(M.Phi.subs(:,1), ...
                                    M.Phi.subs(:,2), ...
                                    M.Phi.subs(:,3), ...
                                    M.Phi.vals, ...
                                    M.DictSig, ...
                                    w_VL, ...
                                    nTheta, ...
                                    nVoxels), ...
                                    size(measured));
rmse_wVL        = sqrt(mean((measured - predicted_VL).^2,1));

% We normalize the rmse by the S0 (non-diffusion measurement)
switch varargin
    case {length(varargin)==1, 'norm'}
      rmse_woVL = rmse_woVL./S0';
      rmse_wVL  = rmse_wVL./S0';
    otherwise
      disp('[feComputeVirtualLesion] Error not normalized by S0')
end

end

% DEBUGGING Plot path, path-neightboord and voxels
% figure('name','path-nehighborhood and tract','color','w')
% hold on
% fg = fe.fg.fibers;
% for ii = 1:length(ind_tract)
%     plot3(fg{ind_tract(ii)}(1,:),fg{ind_tract(ii)}(2,:),fg{ind_tract(ii)}(3,:),'r-');
% end
% view(90,0)
%
% ind_nnz = find(w);
% %tmp_ind = randsample(ind2,200);
% tmp_ind = ind2;
% tmp_ind = intersect(tmp_ind,ind_nnz);
% for ii = 1:length(tmp_ind)
%  plot3(fg{tmp_ind(ii)}(1,:),fg{tmp_ind(ii)}(2,:),fg{tmp_ind(ii)}(3,:),'b-')    
% end
% drawnow;
%
% c = fe.roi.coords - 1.5*ones(size(fe.roi.coords));
% c = c(voxel_ind,:);
% %tmp_c_ind = randsample(1:size(c,1),200);
% tmp_c_ind = 1:size(c,1);
% plot3(c(tmp_c_ind,1),c(tmp_c_ind,2),c(tmp_c_ind,3),'co')

