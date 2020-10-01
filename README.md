### eicsmeardetectors
Collection of smearing scripts for eic-smear.

Legend:
* **Name:** Short descriptor. In general, the corresponding script will start with `Smear` and the detector creation function with `Build`, replacing spaces and periods in the optional version number with underscores. For example, "MatrixDetector 1.2" corresponds to ```SmearMatrixDetector_0_1.cxx```, which implements the function ```BuildMatrixDetector_0_1()```.
* **Min. version:** Recommended or required minimal version of ```eic-smear```. This is usually the eic-smear version at the time the parameterization was added. It is recommended to always use the most current version though; on-going development and bug fixes mean that identical results cannot be guaranteed between different eic-smear versions.
Recent eic-smear versions (since 1.1.0) include an executable that displays the installed version number using
```
eic-smear -v
```
When in doubt, please contact your system administrator.


#### Official parameterizations ####
These are recommended for Yellow Report work. This collection will grow as the detector matrix gets upgraded and concrete designs are parameterized.

<!-- Fill up with  &nbsp; where needed to stop silly line breaks -->

|Name| min. version | Details and Comments |
| --- | --- | --- |
|MatrixDetector 0.1 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | 1.0.4 | Based on the [Detector Matrix](https://physdiv.jlab.org/DetectorMatrix) from June 16 2020. This is the most up-to-date version released by the YR DWG. Note that the HandBook version made assumptions on constant terms and angular resolution where they were/are missing in the matrix. The implementation here explicitly does not do so, and aims to be completely faithful (with perfect angular resolution) to the configuration available online .|
|HandBook 1.2| 1.0.4 | Based on [HANDBOOK v1.2 (Feb 20, 2020)](http://www.eicug.org/web/content/detector-rd) with adaptations where information was incomplete or missing. The MatrixDetector should be considered as the official version. |
|PerfectDetector| 1.0.4 | Perfect detection and PID in \|&eta;\| < 15 |


#### WG additions ####
These are derived from an official detector, customized or extended for specific working group needs.
|Name| min. version | Details and Comments |
| --- | --- | --- |
|MatrixDetector 0.1 with Far Forward detectors | 1.1.0 | Based on the Detector Matrix from June 16 2020 with additional ZDC, B0, and Roman Pots, as found in the [Detector Forward-IR Wiki](https://wiki.bnl.gov/eicug/index.php/Yellow_Report_Detector_Forward-IR). The ZDC only accepts neutrons and photons by default. The Build function accepts the beam momentum per nucleon as an integer parameter. Only 275, 100, 41 (e+P), and 135 (e+D) are accepted. These are ROUGH approximations only!|
|MatrixDetector 0.1 with Barrel TOF | 1.1.1 | Incorporated tofBarrel from https://gitlab.com/preghenella/pid. Based on the Detector Matrix from June 16 2020. This is under active development and not intended to be used widely yet.|

#### Unofficial parameterizations ####
This is a collection of existing parameterizations in various states. They can serve as placeholders and examples until fresh parameterizations are created, approved, and added to the official list.

|Name| min. version | Details and Comments |
|---| ---| --- |
|JLEIC 0.1| 1.1.0 | Adapted from this [stand-alone implementation](https://gitlab.com/eic/escalate/ejana/-/blob/master/src/plugins/reco/eic_smear/SS_DetectorJleic_v1_0_2.cc) in escalate. |
|BeAST 0.1| 1.0.4 | Multi-purpose version of BeAST 0.0, with otherwise the same caveats |
|BeAST 0.0| 1.0.3 | Similar to the version used in [arXiv:1702.00345](https://arxiv.org/abs/1708.01527) and [arXiv:1911.00657](https://arxiv.org/abs/1911.00657). Some shortcuts and tricks were used, making this not suitable for general use. |
|ePHENIX 0.0| 1.0.3 | An example implementation from 2014. Note that it requires compilation; instructions are in the script. Some details are in the comments at the bottom of the file.|
|eSTAR 0.0| 1.0.3 | An example implementation from 2014.  Based on parameterizations given in [Zhangbu Xu's talk (slide 5)](https://wiki.bnl.gov/conferences/index.php/January_2014) |
|ZEUS 0.0| 1.0.3 | An example implementation from 2014. Based on JHEP05 (2009) 108|
|STAR 0.0| 1.0.3 | Some details are in the comments at the bottom of the file. |

#### Basic interactive usage ####

Transformation and smearing can be performed interactively in ROOT as follows:

##### Generate EicTree #####

```
root [] gSystem->Load("libeicsmear");
root [] BuildTree ("ep_hiQ2.20x250.small.txt.gz",".",-1, "");
Processed 10000 events containing 346937 particles in 6.20576 seconds (0.000620576 sec/event)
```
BuildTree accepts the name of an input file, the output directory, and
the number of events to generate (-1 for all). The final string argument
optionally accepts the name of a logfile created with the MC generator by using the ">" redirection. If it is not provided, the macro will look for a file with the same name as the input file but with the ending ".log". It will search in the same directory or, if the txt file is in a directory named "example/TXTFILES", in a directory named "example/LOGFILES".
These log files contain additional information on generated cross section, events, and for some generators trials which are needed to calculate the total cross section.

IMPORTANT: The file type is by default assumed to be pythia6. For
other files, please make sure to include the generator name in the
filename. Currently accepted are hepmc as well as pythia, pepsi, lepto, rapgap, djangoh, beagle, milou, sartre, simple. If the filename ends with ".gz" or ".zip", it will be decompressed on the fly, so there is no need to unzip it.


#### Smear the tree
```
root [] gSystem->Load("libeicsmear")
root [] .L SmearMatrixDetector_0_1.cxx // Assuming you copied this here
root [] SmearTree(BuildMatrixDetector_0_1(),"ep_hiQ2.20x250.small.root", "smeared.root",-1)
/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/
/  Commencing Smearing of 10000 events.
/-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-/
```
SmearTree accepts the name of an input file, the output file, and
the number of events to smear (-1 for all);

* IMPORTANT: The file type is by default assumed to be pythia6. For
other files, please make sure to include the generator name in the filename.

* You may note many warnings of the form
```
Warning in Kinematics, bounded(): x (or y) = 1.02488 is outside [0,1]
To disable this warning, set erhic::DisKinematics::BoundaryWarning=false;
```
These are largely benign - a smeared electron or missing/smeared hadronic information can lead to "unphysical" x values. The warning is easily suppressed using
```
root [] erhic::DisKinematics::BoundaryWarning=false;
root [] SmearTree(BuildMatrixDetector_0_1(),"ep_hiQ2.20x250.small.root", "smeared.root",-1)
```
but such warnings can serve as canaries in a coal mine while experimenting with a smearing script or an MC generator, especially if "nan" or "inf" values are produced, and are therefore turned on by default.

##### Analyze the Tree #####

(Suppressing the root prompts for easier copy/paste):
```
root -l
gSystem->Load("libeicsmear");
TFile mcf ("ep_hiQ2.20x250.small.root"); // truth
TTree* mc=(TTree*)mcf.Get("EICTree");
mc->AddFriend("Smeared","smeared.root"); // befriend

// Setup Input Event Buffer
erhic::EventMC* inEvent(NULL);
Smear::Event* inEventS(NULL);
mc->SetBranchAddress("event",&inEvent);
mc->SetBranchAddress("eventS",&inEventS);

// histo and event loop
TH2D* EEprime = new TH2D("EEprime", ";E;Eprime", 100, 0, 20, 100, 0, 20);
for(long iEvent=0; iEvent<mc->GetEntries(); iEvent++){    
    //Read next Event
    if(mc->GetEntry(iEvent) <=0) break;
    // Loop over Particles
    for(int j=0; j<inEventS->GetNTracks(); j++){
      if ( j<3 ) continue;       // Skip beam
      const Smear::ParticleMCS* inParticleS = inEventS->GetTrack(j); // Smeared Particle      
      const Particle* inParticle = inEvent->GetTrack(j); // Unsmeared Particle
      if(inParticleS == NULL) continue;
      EEprime->Fill(inParticle->GetE(), inParticleS->GetE());
   }
}
EEprime->Draw("colz");
```

#### Compilation, Tests and Examples ####

To create a library and tests and examples follow these instructions:
```
mkdir build
cd build
cmake -DBUILD_TESTS=ON -DCMAKE_INSTALL_PREFIX=$EICDIRECTORY ..
make
```
The ```CMAKE_INSTALL_PREFIX``` is optional and only needed if you would like to install the library using
```
make install
```

Depending on your environment, you may need to give cmake a hint where to find the eic-smear configuration. If it has been installed into $EICDIRECTORY, do this (adapt as appropriate for other locations, and keep in mind ```setenv A B``` in csh corresponds to ```export A=B``` in bash):
```
setenv CMAKE_PREFIX_PATH $EICDIRECTORY/cmake:$CMAKE_PREFIX_PATH
```

To properly use the library and override an existing system-wide one, prepend the location to adapt the LD_LIBRARY_PATH
```
setenv LD_LIBRARY_PATH $EICDIRECTORY:$LD_LIBRARY_PATH
```
or
```
setenv LD_LIBRARY_PATH .:$LD_LIBRARY_PATH
```

Also, in case you see this error:  'eicsmear/smear/Detector.h' file not found, make sure that you set the include path like so, before compiling:

```
export ROOT_INCLUDE_PATH=path-to-eic-smear-install/include
```
A variety of tests are generated from the tests directory:

```
./test_simple_buildtree
```
will read a (provided) e+D BeAGLE file.

```
./particlegun
```
is a customizable particle gun that creates a few simple histograms
and plots to see and test the acceptance dependence of smearing.

And of course, you can now use all detector scripts from the ROOT interface with
one line as well.
```
root [] gSystem->Load("libeicsmear")
root [] gSystem->Load("libeicsmeardetectors")
root [] SmearTree(BuildByName("MATRIX"),"ep_hiQ2.20x250.small.root")
```

A wrapper in eic-smear allows to start ROOT with the libraries loaded and displays
version information as well as the library locations.
```
$ eic-smear
Using eic-smear version: 1.1.0-rc1
Using these eic-smear libraries :
/Users/kkauder/software/lib/libeicsmear.dylib
/Users/kkauder/software/lib/libeicsmeardetectors.dylib
eic-smear [0] SmearTree(BuildByName("MATRIX"),"ep_hiQ2.20x250.small.root")
```

It can also be used for simple one liners:
```
echo 'BuildTree ("ep_hiQ2.20x250.small.txt.gz");SmearTree(BuildByName("MATRIX"),"ep_hiQ2.20x250.small.root")' | eic-smear
```

#### A canonic example ####

When tests are built, a particularly useful example is```qaplots```.
This starts from (provided) examples of text MC output, builds the
EicTree, smears it, and runs a simple analysis on the result. The
current set of QA plots will grow, but the source file
```
tests/qaplot.cxx
```
and associated cmake configuration in ```CMakeLists.txt```
is intended to be readable and give insight how to generate compiled
code.

You can use them with the supplied test files which cmake copies into
your build directory:
```
 ./qaplots -i ep_lowQ2.20x250.small.txt.gz -det matrix
```
or  
```
./qaplots -i ep_hiQ2.20x250.small.txt.gz -det matrix
```

The first set of event-wise observables (y, x, and Q^2 using
three different methods) is set up to loosely compare to the plots on
p. 88f in the eRHIC design study,  http://arxiv.org/pdf/1409.1633.pdf,
but be aware of statistics limitations and specific settings (e.g.,
20x250 GeV e+P) in the test files.

You may wish to combine the two test files. In order to do so, the
repeated header in the second (or any other subsequent) file needs to
be removed. This can be achieved like this:
```
cp ep_lowQ2.20x250.small.txt.gz combo.txt.gz
gunzip combo.txt.gz
gunzip ep_hiQ2.20x250.small.txt
tail -n +7 ep_hiQ2.20x250.small.txt >> combo.txt
gzip combo.txt
gzip ep_hiQ2.20x250.small.txt
```

The result of
```
./qaplots -i combo.txt.gz -det handbook
```
should then reproduce (within the limits of architecture-dependent
random number details) the provided small reference plots in ```epref_qaplotsHANDBOOK.pdf```.


#### Anatomy of a Smearer ####

A "detector" is constructed as follows. For additional details,
please also see the provided source files and comments therein.

```
// ... omitted some includes and helpers
Smear::Detector BuildMyDetector() {

// Create the detector object to hold all devices
   Smear::Detector det;
  // The detector will calculate event kinematics from smeared values
  det.SetEventKinematicsCalculator("NM JB DA");
```

* Set up a device that smears. In this case, momentum is smeared
  * in eta = -3.5 --  -2.5,
  * with sigma_p/p = 0.1 % p+ 2.0 %,
  * accepting all charged particles.

```
  // Tracking
  // eta = -3.5 --  -2.5
  // sigma_p/p ~ 0.1% p+2.0%
  Smear::Acceptance::Zone TrackBack1Zone(ThetaFromEta ( -2.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device TrackBack1P(Smear::kP, "sqrt( pow ( 0.001*P*P, 2) + pow ( 0.02*P, 2) )");
  TrackBack1P.Accept.AddZone(TrackBack1Zone);
  TrackBack1P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBack1P);
  ```

* IMPORTANT: For more realistic representation of electrons, it may
  make sense to only accept hadrons and create a separate device for
  electrons that represents the combined information from multiple
  detectors, e.g.:

```
  TrackBack1P.Accept.SetCharge(Smear::kCharged);
  TrackBack1P.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(TrackBack1P);

  Smear::Device EMTrackBack1P(Smear::kP, "sqrt(  0.01 * pow(p.Pt(), 2) + 0.005 * p.Pt())");
  EMTrackBack1P.Accept.AddZone(TrackBack1Zone);
  EMTrackBack1P.Accept.SetCharge(Smear::kCharged);
  TrackBack1P.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EMTrackBack1P);
```

* Continue on, adding $`\phi`$ and $`\theta`$ devices, calorimetry,
etc.

* In case you want to avoid consistency checks (not recommended), you can
add

```
det.SetLegacyMode(true); // Debug and experimentation only
```

* Finally, after adding all desired devices return the complete
detector.

```
  return det;
}
```


* Formulas are based on ROOT::TFormula and accept P, Phi, Theta,
E, Pt and Pz. The logic also supports Px and Py, but this
isn't included in the parser yet; please contact the authors if you
need this functionality added.

##### IMPORTANT NOTES: #####

* If you want to have an unsmeared value in the smeared tree, use a
perfect device, e.g:
```
Smear::Device TrackTheta(Smear::kTheta, "0");
```

* Once any variable is smeared, all remaining fields are initialized
to 0, meaning eic-smear treats them as measured with a zero value.
In recent versions of eic-smear, one can determine whether a property is smeared (and just happens to be 0) via flags of the form
```
if ( inParticleS->IsThetaSmeared() ) { \\ physical, treat as measured };
```
For smeared trees created with older versions, it is up to the user to catch this behavior with lines of
the form
```
if ( fabs ( inParticleS->GetTheta()) > 1e-8 ) { \\ physical, treat as measured }
```
See also some more details in SmearHandBook_1_2.cxx
regarding calculation of smeared kinematic variables.

* Due to the limitations of TFormula only four different
  variables can be used at a time (because internally they get
  translated dynamically into the four available free variables
  x,y,z,t).

* If two devices of the same type have overlapping acceptance, an error will be generated at runtime if a particle is smeared by both. There is no realistic way to automatically detect colliding acceptances for arbitrary acceptance zones.

* Recent changes make it a requirement to always have "location" information. In practice, this means:
  * provide exactly 0 or exactly three smeared momentum components. The others will then be automatically constructed. If phi or theta resolution is unknown, this can be of the form ```Smear::Device phi(Smear::kPhi,"0");```
  * In the case of no momentum smearing (i.e., only energy is smeared), angular smearers still need to be provided, now representing phi and theta as determined by the calorimeter.
  * To support legacy smearing scripts, this behavior can be turned off via ```Detector::SetLegacyMode(true);```
  Note that in this case Pt and Pz smearing will no longer work (as that was the status when the legacy scripts were created).

### More Examples ###
Additional usage examples can be found at https://github.com/eic/eicsmear-jetexample. They are meant to help in any analysis, not just jet-based ones.
