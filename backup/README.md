This directory conains methods for implementing various beam effects which 
will be present at the EIC (crossing angle, divergence, crabbing kicks, etc) 
in Pythia8 using the native BeamShape functionality. 

eicBeamShape.h/cxx: Definition and implementation of beam effects

PythiaBeamShape.cxx: Main code to run Pythia8 with beam effects implemented and generate various diagnostic plots in ROOT. Contains several kinematic cuts appropriate for jet / SIDIS type analyses.

PythiaBeamShapeHepMC.cxx: Code to run Pythia8 and pass output as HepMC file. No analysis is performed and no cuts are introduced.

dis.steer: Dummy steering file to set Pythia8 configuration when running tests

steerFiles/: Directory containing steering files for the beam energy combinations and configurations (high divergence or high acceptance) listed in the CDR. NC DIS events are simulated with Q2 > 10 GeV^2


To run PythiaBeamShape.cxx:

make

./runBeamShape.cxx steerFile configFlag hadEnergy lepEnergy x-Angle out.hist.root

Note: x-Angle (crossing angle) needs to be input in radians, ie 0.025 for a 25 milliradian crossing angle. configFlag = 1 is for hiDiv, configFlag = 2 is for hiAcc, and configFlag = 3 is for when running eA energies (18x110, 10x110, 5x110, and 5x41).

Note: At IP6, we have instituted a right-handed coordinate system with the x-direction pointing to the ring center, the y-direction pointing to the sky, and the z-direction pointing in the direction of the hadron beam. As the beams enter IP6 from ring inside, the hadron beam points to the negative x positive y direction. So for IP6, one must input the crossing angle with a negative sign: x-Angle = -0.025.