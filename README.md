# active-information-engine

This repository contains the code, data, figures and figure scripts that generate supporting numerical results reported in:

Garcia-Millan*, Schüttler*, Cates, and Loos, Optimal closed-loop control of active particles and a minimal information engine, arXiv:2407.18542 [cond-mat.stat-mech] (2024)

Schüttler, Garcia-Millan, Cates, and Loos, Active particles in moving traps: minimum work protocols and information efficiency of work extraction, arXiv:2501.18613 [cond-mat.stat-mech] (2025)

To run the code, copy and paste the following compilation and execution commands:

cc -O3 -Wall -o eng rnt_engine.c -lgsl

GSL_RNG_SEED=$RANDOM ./eng

cc -O3 -Wall -o eng rnt_engine_distrW.c -lgsl

GSL_RNG_SEED=$RANDOM ./eng

If you use this software or data, please cite it as:

@unpublished{Garcia-MillanSchuettlerETAL:2024,
  title={Optimal closed-loop control of active particles and a minimal information engine},
  author={Garcia-Millan, Rosalba and Sch{\"u}ttler, Janik and Cates, Michael E and Loos, Sarah AM},
  eprint={arXiv:2407.18542},
  url={ https://arxiv.org/pdf/2407.18542 },
  year={2024}
}

@unpublished{SchuettlerETAL:2025,
  title={Active particles in moving traps: minimum work protocols and information efficiency of work extraction},
  author={Sch{\"u}ttler, Janik and Garcia-Millan, Rosalba and Cates, Michael E and Loos, Sarah AM},
  eprint={arXiv:2501.18613},
  url={ https://arxiv.org/pdf/2501.18613 },
  year={2025}
}
