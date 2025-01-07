classdef userSettings
    %holds the changable parameters for homogenization + nozzle
    %optimization 

    properties
        %% mesh properties/dimensions        
        nelx = 100;
        %nely and nelz are calculated in top3d_anisotropy using
        %getDims(obj)
        tFlange = 0.5;
        tWeb= 0.9;
        volfrac = 0.5;
        isPassive = 1; %==1 if passive region is turned on
        %% Filtering
        beta = 20; %initial
        beta_max = 20; %continuation isn't being used right now. Check line 197 if you want to turn it on
        d_beta = 1.1; %change in beta

        eta = 3; %initial
        eta_max = 10;
        d_eta = 1;
        rho_min = 1e-9;
        %% Filter Geometry
        %dim for circle filter
        rmin = 2;
        %% material properties:
        Emin = 1e-9;
        nu = 0.3;
        E1 = 10000;
        E2 = 1000;
        E3 = 1000;

        %%
        isotropy = 1; %==0 if isotropic, ==1 if anisotropic

        maxloop = 10000;    % Maximum number of iterations
        tolx = 0.0001;      % Terminarion criterion
        displayflag = 0;  % Display structure flag. if ==1 plots every 10 iterations
    end

    methods
        function [nely, nelz] = getDims(obj)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            nely = round(obj.nelx/10);
            nelz = round(nely/2);
        end

        function [filename] = getFilename(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if obj.isotropy == 1
                iso = "_ani_";
            else
                iso = "_iso_"
            end
            vol = "_" + string(obj.volfrac) + "";

            filename = string(obj.nelx) + vol + iso + 'beam.mat';
        end
    end
end