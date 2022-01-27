function plot_network(trails, phi, grids, params, repeat_index, field_skip, field_sf)
%% Plots the network and the field after repeat_index>0 number of repeats.

if nargin < 6
    field_skip = 3;
end
if nargin < 7
    field_sf = 2.5;
end

figure
hold on
% Draw the vasculature.
toPlot = double(trails{repeat_index}~=0); toPlot = full(toPlot); toPlot(~toPlot | ~grids.within_bounds_mask) = NaN;
pcolor(grids.mesh_X,grids.mesh_Y,toPlot)
trailColor = [0,0,0];
colormap([trailColor])
shading flat
x = grids.phi_mesh_X(:); y = grids.phi_mesh_Y(:);
basis = basis_func(x,y,params);
lin_inds = 1 : numel(phi.phi_0);
phiNow = reshape(eval_phi_basis(lin_inds, phi.phi_0(lin_inds(:)), basis, phi.weights{repeat_index}),size(phi.phi_0));
dx = grids.phi_grid_X(2) - grids.phi_grid_X(1); dy = grids.phi_grid_Y(2) - grids.phi_grid_Y(1);
q = quiver([grids.phi_mesh_X(1:field_skip:end,1:field_skip:end);grids.phi_mesh_X(1:field_skip:end,1:field_skip:end)], [grids.phi_mesh_Y(1:field_skip:end,1:field_skip:end);grids.phi_mesh_Y(1:field_skip:end,1:field_skip:end)], field_sf*dx*[cos(phiNow(1:field_skip:end,1:field_skip:end));-cos(phiNow(1:field_skip:end,1:field_skip:end))]/2, field_sf*dy*[sin(phiNow(1:field_skip:end,1:field_skip:end));-sin(phiNow(1:field_skip:end,1:field_skip:end))]/2,...
    0,'ShowArrowHead','off','LineWidth',1.5,'Color',0.6*[1,1,1]);
axis equal
xlim(params.outerRadius*[-1,1])
ylim(params.outerRadius*[-1,1])
axis off
% Draw the annulus.
theta = linspace(0,2*pi,1e3);
plot(params.innerRadius*cos(theta), params.innerRadius*sin(theta),'LineWidth',2,'Color','black')
plot(params.outerRadius*cos(theta), params.outerRadius*sin(theta),'LineWidth',2,'Color','black')