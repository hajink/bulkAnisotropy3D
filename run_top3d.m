function run_top3d()
%simple wrapper to run top3d by Liu and Tovar
var = userSettings(); %object class that hold all of the user-defined parameters

%can iterate through runs like the following: 
%changes volume and isotropy
vols = [0.355];
isos = [0 ];

for i = 1:length(vols)
    for j = 1:length(isos)
        var.volfrac = vols(i);
        var.isotropy = isos(j);
        top3d_anisotropy(var)
    end
end


return