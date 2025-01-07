A neat(er) and shareable code to generate 3d topology optimized structures with anisotropic base material assumptions[1]. Currently setup to generate I-beam-like designs for 3-point bending.\

The base of the code has been updated to be based on https://www.top3d.app/ [2]. The largest modification is to the generation of the stiffness matrix KE. New functions are in the folder titled *stiffnessGeneration*. The function that plays the largest role is *D_3d.m*, which initializes and populates the constitutive matrix. MMA[3] is implemented as the optimizer. 

To run, user may simply run run_top3d(). To alter user-input variables/parameters, modifications can be made to the *userSettings* object declared in line 3 of run_top3d(). Users are encouraged to read the class definition to see what variables and subfunctions are modifiable. 

References: \
[1] Kim, H., and Carstensen, J.V. Experimental Investigation of Topology-Optimized Beams with Isotropic and Anisotropic Base Material Assumptions, ASME IDETC 2022, St. Louis (MO) pp. 1-11, August 2022. https://asmedigitalcollection.asme.org/IDETC-CIE/proceedings/IDETC-CIE2022/86236/V03BT03A033/1150428 \
[2] Liu, K., Tovar, A. An efficient 3D topology optimization code written in Matlab. Struct Multidisc Optim 50, 1175–1196 (2014). https://doi.org/10.1007/s00158-014-1107-x \
[3] Svanberg K (1987) The method of moving asymptotes-a new method for structural optimzation. Int J Numer Methods Eng 24:359–373. https://doi.org/10.1002/nme.1620240207
