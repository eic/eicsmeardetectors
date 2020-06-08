### eicsmeardetectors
Collection of smearing scripts for eic-smear.

Legend:
* **Name:** Short descriptor. In general, the corresponding script will start with `Smear` and the detector creation function with `Build`, replacing spaces and periods in the optional version number with underscores. For example, "Handbook 1.2" corresponds to ```SmearHandBook_1_2.cxx```, which implements the function ```BuildHandBook_1_2()```.
* **Min. version:** Recommended or required minimal version of ```eic-smear```. This is usually the eic-smear version at the time the parameterization was added. It is recommended to always use the most current version though; on-going development and bug fixes mean that identical results cannot be guaranteed between different eic-smear versions.
Future eic-smear versions will include an executable that displays the installed version number using
```
eic-smear -v
```
When in doubt, please contact your system administrator.


#### Official parameterizations ####
These are recommended for Yellow Report work. This collection will grow as the detector matrix gets upgraded and concrete designs are parameterized.

XXX

|Name| min. version | Details and Comments |
|---| ---| --- |
|MatrixDetector 0.1| 1.0.4 | Based on the [Detector Matrix](https://physdiv.jlab.org/DetectorMatrix) from June 8 2020|
|HandBook 1.2| 1.0.4 | Based on [HANDBOOK v1.2 (Feb 20, 2020)](http://www.eicug.org/web/content/detector-rd) with adaptations where information was incomplete or missing|
|PerfectDetector| 1.0.4 | Perfect detection and PID in \|&eta;\| < 15 |

#### Unofficial parameterizations ####
This is a collection of existing parameterizations in various states. They can serve as placeholders and examples until fresh parameterizations are created, approved, and added to the official list.

|Name| min. version | Details and Comments |
|---| ---| --- |
|BeAST 0.1| 1.0.4 | Multi-purpose version of BeAST 0.0, with otherwise the same caveats |
|BeAST 0.0| 1.0.3 | Similar to the version used in [arXiv:1702.00345](https://arxiv.org/abs/1708.01527) and [arXiv:1911.00657](https://arxiv.org/abs/1911.00657). Some shortcuts and tricks were used, making this not suitable for general use. |
|ePHENIX 0.0| 1.0.3 | An example implementation from 2014. Note that it requires compilation; instructions are in the script. Some details are in the comments at the bottom of the file.|
|eSTAR 0.0| 1.0.3 | An example implementation from 2014.  Based on parameterizations given in [Zhangbu Xu's talk (slide 5)](https://wiki.bnl.gov/conferences/index.php/January_2014) |
|ZEUS 0.0| 1.0.3 | An example implementation from 2014. Based on JHEP05 (2009) 108|
|STAR 0.0| 1.0.3 | Some details are in the comments at the bottom of the file. |
