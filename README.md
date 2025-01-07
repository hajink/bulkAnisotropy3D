The base of the code is from https://www.top3d.app/ [1]. The largest modification is to the generation of the stiffness matrix KE. New functions are in the folder titled *stiffnessGeneration*. The function that plays the largest role is *D_3d.m*, which initializes and populates the constitutive matrix. MMA[2] is implemented as the optimizer. 

To run, user may simply run run_top3d(). To alter user-input variables/parameters, modifications can be made to the *userSettings* object declared in line 3 of run_top3d(). Users are encouraged to read the class definition to see what variables and subfunctions are modifiable. 

References: \
[1]Liu, K., Tovar, A. An efficient 3D topology optimization code written in Matlab. Struct Multidisc Optim 50, 1175–1196 (2014). https://doi.org/10.1007/s00158-014-1107-x \
[2] Svanberg K (1987) The method of moving asymptotes-a new method for structural optimzation. Int J Numer Methods Eng 24:359–373. https://doi.org/10.1002/nme.1620240207
