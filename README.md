# active-information-engine

This repository contains the code, data, figures and figure scripts that generate supporting numerical results reported in:

Garcia-Millan*, Schüttler*, Cates, and Loos, Optimal closed-loop control of active particles and a minimal information engine, Phys. Rev. Lett. 135, 088301 (2025)

https://doi.org/10.1103/fbgp-qpvv

[preprint arXiv:2407.18542]


Schüttler, Garcia-Millan, Cates, and Loos, Active particles in moving traps: Minimum work protocols and information efficiency of work extraction, Phys. Rev. E 112, 024119 (2025)

https://doi.org/10.1103/4q4f-1dpx

[preprint arXiv:2501.18613]


To run the code, copy and paste the following compilation and execution commands:

cc -O3 -Wall -o eng rnt_engine.c -lgsl

GSL_RNG_SEED=$RANDOM ./eng

cc -O3 -Wall -o eng rnt_engine_distrW.c -lgsl

GSL_RNG_SEED=$RANDOM ./eng

If you use this software or data, please cite it as:

@article{Garcia-MillanSchuettlerETAL:2025,

  title={Optimal closed-loop control of active particles and a minimal information engine},
  
  author={Garcia-Millan, Rosalba and Sch{\"u}ttler, Janik and Cates, Michael E and Loos, Sarah AM},
  
  journal={Phys. Rev. Lett.},
  
  volume={135},
  
  number={8},
  
  pages={088301},
  
  year={2025},
  
  publisher={APS}

}


@article{SchuettlerETAL:2025,

  title={Active particles in moving traps: Minimum work protocols and information efficiency of work extraction},
  
  author={Sch{\"u}ttler, Janik and Garcia-Millan, Rosalba and Cates, Michael E and Loos, Sarah AM},
  
  journal={Phys. Rev. E},
  
  volume={112},
  
  number={2},
  
  pages={024119},
  
  year={2025},
  
  publisher={APS}

}
