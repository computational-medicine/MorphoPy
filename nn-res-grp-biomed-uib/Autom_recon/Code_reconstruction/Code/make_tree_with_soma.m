function make_tree_with_soma(somaContour,tree,swcfile)

% generate tree toolbox structure for export
Npoints = length(somaContour.X);
somacontourtree.dA = sparse(Npoints,Npoints);
somacontourtree.X = somaContour.X;
somacontourtree.Y = somaContour.Y;
somacontourtree.Z = somaContour.Z;

somacontourtree.D = ones(Npoints,1)*0.1; % diameter of contour set to 0.1. Only affects the line size of the contour on screen.
% set type to 1=soma
somacontourtree.R = ones(Npoints,1);
somacontourtree.rnames{1} = 'soma';
somacontourtree.rnames{2} = '';
somacontourtree.rnames{3} = 'dendrite';
somacontourtree.dA = sparse(Npoints,Npoints);
for n=2:Npoints
   somacontourtree.dA(n,n-1) = 1;
end
somacontourtree.name = 'somacontour';

% tree should be a single origin point, with all subtrees attached to it.
% so attach individual subtrees to soma contour one by one (last point of
% soma contour)
subtree_ind = find(tree.dA(:,1))';
treewithSoma = somacontourtree;
for n = subtree_ind
[~, subtree] = sub_tree(tree, n);
subtree.R = subtree.R*0+3;    % set all subtree region codes to 3 ('dendrites')
subtree = rmfield(subtree,'rnames'); % remove region names (strictly speaking not necessary here)
% attach to last point of soma contour (Npoints) to keep the contour in one place in the exported file
% use '-r' flag to omit any updates of area codes, which messes up the codes when exporting
treewithSoma = cat_tree(treewithSoma,subtree,Npoints,[],'-r');   
%treewithSoma.R = treewithSoma.R';
end

swc_tree(treewithSoma,swcfile);

end