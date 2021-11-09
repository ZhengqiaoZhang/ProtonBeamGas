To run BeamGasEvent.cxx:

make

./runBeamGasHepMC.exe steerFiles/dis_eicBeam_hiAcc_18x275 2 275 18 -0.025 out.hist.root out.hepmc

Note: x-Angle (crossing angle) needs to be input in radians, ie 0.025 for a 25 milliradian crossing angle. configFlag = 1 is for hiDiv, configFlag = 2 is for hiAcc, and configFlag = 3 is for when running eA energies (18x110, 10x110, 5x110, and 5x41).

Note: At IP6, we have instituted a right-handed coordinate system with the x-direction pointing to the ring center, the y-direction pointing to the sky, and the z-direction pointing in the direction of the hadron beam. As the beams enter IP6 from ring inside, the hadron beam points to the negative x positive y direction. So for IP6, one must input the crossing angle with a negative sign: x-Angle = -0.025.
