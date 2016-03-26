#define GRID 0

#if GRID == 1
#define SOLVER_GRID_UG

#else
#define SOLVER_GRID_NON_UG

#endif

#define SOLVER 1

#if SOLVER == 1
#define POISSON_SOLVER_JACOBI
#endif

#if SOLVER == 2
#define POISSON_SOLVER_GS
#endif

#if SOLVER == 3
#define POISSON_SOLVER_GSOR
#endif

#define MAX_STRING_LENGTH 80
