#define GRID 2

#if GRID == 1
#define SOLVER_GRID_UG

#else
#define SOLVER_GRID_NON_UG

#endif

#define SOLVER 2

#if SOLVER == 1
#define POISSON_SOLVER_JACOBI
#endif

#if SOLVER == 2
#define POISSON_SOLVER_GS
#endif

#if SOLVER == 3
#define POISSON_SOLVER_GSOR
#define omega 1.1
#endif

#define MAX_STRING_LENGTH 80


#define Nxb 24
#define Nyb 24


#define D_xmin -0.5
#define D_ymin -0.5


#define D_xmax 0.5
#define D_ymax 0.5


#define HK 2
#define HD 2
 

#define MaxIt 1500


#define PRES_VAR 1
#define CENT_VAR 1

#define VELC_VAR 1
#define FACE_VAR 1

#define inRe 0.001

